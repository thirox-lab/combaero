# Extraction: Lau, Kukreja and McMillin (1990), V-shaped rib arrays

**Status: CONFIRMED for the 90 deg baseline (case 1a/-90), which is the only
configuration an implemented set can be compared against. Table 2 is confirmed
cell by cell on three independent channels -- all 42 numbers, no exceptions.
The V-shaped cases are transcribed and, having now been verified to the same
standard, are held back only on the modelling question in L1, not on any doubt
about the numbers.**

Tracked by #385, under #333.

> Lau, S.C., Kukreja, R.T. and McMillin, R.D. (1990). "Effects of V-shaped rib
> arrays on turbulent heat transfer and friction of fully developed flow in a
> square channel." *Int. J. Heat Mass Transfer* 34(7), 1605-1616.
> `sources/heat_transfer/rib/ribs_lau_1990.pdf` (gitignored, copyrighted).

**The first rib source outside Han**, and a primary paper rather than a
textbook reprint. #333 named this gap explicitly: agreement with the curve a
correlation was fitted to is a consistency check, not a validation.

It also arrives **tabulated**. Table 2 prints the roughness functions as
closed forms, so nothing about the correlation needed digitising -- which is
the preference `extractions/README.md` states.

---

## For the reviewer

| # | what needs judgement | what turns on it |
|---|---|---|
| **L1** | **Is a Lau V-rib `RibCorrelationSet` wanted?** Table 2 gives `R`/`G`/`G_bar,Lau` coefficients for V-45, V-60, V-120 and V-135, a shape `han_park_1988_angled` cannot express (it carries angle only, which is why Han figure 4.51's own V classes are deliberately unscored). | A code task shaped like #334, not a data task. **The coefficients are now three-channel verified**, so this is purely a modelling decision -- nothing further needs reading. |
| **L2** | **Is the ~13% `R` gap acceptable as a recorded disagreement?** **Now SCORED rather than narrated** (#339 plan 1a): MAE 11.3%, bias -11.3%, and **0 of 9 points inside Lau's own stated +/-10.9% friction band**. The harness reports it on every run instead of a doc arguing it. Definitional causes were checked and ruled out (below). Two labs differ by ~12% in ribbed-wall friction for nominally the same 90 deg configuration. | Nothing breaks either way, but it bounds how much trust the friction side of `han_1988_orthogonal` should carry. A third source would settle whether Han or Lau is the outlier. |

### What I still need sent

Nothing. The scan's right margin is clipped on page 2, so the rib-geometry
sentence was read in fragments -- but the numbers are confirmed by an exact
arithmetic identity (below) and by the reviewer independently.

---

## Geometry and conditions

| quantity | value | |
|---|---|---|
| channel | square, aluminium, 1.52 m long, 7.62 x 7.62 cm | confirmed; reviewer read the same sentence independently |
| `D` | 0.0762 m | square, so `D_h` = side |
| `L/D` | 19.9 | |
| rib height `e` | 4.76 mm | read from a clipped column, confirmed by the identity below |
| **`e/D`** | **0.0625** | |
| ribbed walls | two opposite; two smooth | |
| `p/e` | 10 and 20 | |
| `alpha` | 45, 60, 90, 120, 135 deg | V-shaped; case 1a is 90 deg full ribs |
| `Re_D` | 10 000 - 60 000 | |
| measurement window | `x/D` = 7.66 - 15.16, thermally fully developed | |
| properties | at average bulk temperature | |

**A check that could have failed.** `e/D` was read off a page whose right
margin the scan cuts. It is exact by design, not by rounding: 7.62 cm is
exactly 3 inches and 4.76 mm exactly 3/16 inch, so
`e/D = (3/16)/3 = 1/16 = 0.0625`. An imperial design reproducing an exact
binary fraction is not something a misread digit produces.

**Stated experimental uncertainty** (Kline and McClintock, ref [6]):
`Re` +/-2.9%, Stanton numbers +/-5.8%, friction factor +/-10.9%. This is a
property of the measurement, so it is what populates `uncertainty` in the
committed metadata -- the distinction #389 exists to enforce.

---

## Extracted items

| # | item | content | status |
|---|---|---|---|
| 1 | Eq. (5), roughness Reynolds number | `e+ = (e/D) Re_D [(2 fbar - f_ss)/2]^0.5` | confirmed |
| 2 | Eq. (6), friction roughness function | `R = [(2 fbar - f_ss)/2]^-0.5 + 2.5 ln[2(e/D)] + 2.5` | confirmed |
| 3 | Eq. (7), heat transfer roughness function | `G = [(2 fbar - f_ss)/2]^0.5 / St_r + 2.5 ln[2(e/D)] + 2.5` | confirmed |
| 4 | Eq. (8), average heat transfer roughness function | `G_bar,Lau = [(2 fbar - f_ss)/2]^0.5 / St_avg + 2.5 ln[2(e/D)] + 2.5`. Written `G_bar,Lau` throughout, never bare `G_bar`, so the distinction is visible wherever it is used | confirmed -- **and not Han's `G_bar`, see item 8** |
| 5 | `(2 fbar - f_ss)/2` | the ribbed-wall friction extracted from the channel average, because only two of four walls are ribbed. `fbar` alone is the channel average | confirmed |
| 6 | `St_r` basis | ribbed-wall heat flux divided by the **projected** area, "not including the increased rib surface area" | confirmed |
| 7 | Table 2 | seven cases, `R`/`G`/`G_bar,Lau` each as `a(e+)^b` -- reproduced below | **confirmed cell by cell on three independent channels, all 42 numbers.** Visual read at 300 dpi, macOS Vision OCR at 600 dpi, and the reviewer. The single OCR disagreement was escalated to a high-resolution crop, resolved, and subsequently confirmed independently |
| 8 | **`Gbar` is a four-wall average** | uses `St_avg` over two ribbed and two smooth walls. **So is Han's `G_bar`** -- `han_ribbed.md` item 10's current resolution is `G_bar = 1.2 G`, printed on figure 4.47, with the Prandtl reading recorded there as withdrawn. An earlier version of this document quoted the withdrawn form and inverted the conclusion; see #392 and the section below | confirmed -- the two are the **same quantity**, and Lau's is committed |
| 9 | geometry group | Lau's `2(e/D)`; combaero's `(2 e/D)(2W/(W+H))`. For a square channel `2W/(W+H) = 1`, so they coincide | confirmed by derivation |

### Table 2, as printed

| case | `R` a | `R` b | `G` a | `G` b | `Gbar` a | `Gbar` b |
|---|---|---|---|---|---|---|
| 1a/-90 | 3.674 | -0.0033 | 4.218 | 0.257 | 5.450 | 0.250 |
| 4/-V-45 | 1.605 | 0.0015 | 1.819 | 0.355 | 2.656 | 0.335 |
| 5/-V-60 | 1.232 | 0.0454 | 1.299 | 0.399 | 2.163 | 0.360 |
| 5a/-V-60 | 3.687 | -0.0571 | 1.685 | 0.376 | 2.719 | 0.336 |
| 6/-V-120 | 1.537 | -0.0228 | 1.983 | 0.352 | 3.460 | 0.315 |
| 6a/-V-120 | 3.617 | -0.0734 | 2.739 | 0.324 | 4.351 | 0.291 |
| 7/-V-135 | 2.070 | -0.0384 | 1.992 | 0.362 | 3.685 | 0.313 |

---

## `Gbar`: a hazard that turned out to be the opposite

This section previously argued that Lau's `Gbar` and Han's `G_bar` were
different quantities sharing a symbol, and that committing Lau's would read
as a false confirmation. **That was wrong, and the way it was wrong is worth
keeping.**

It rested on `han_ribbed.md` item 10 in its **withdrawn** form -- the reading
that Han's `G_bar` is the Prandtl-normalised `G Pr^-0.57` (factor 1.2162).
Item 10 had already corrected itself three hundred lines further down: figure
4.47 **prints** `G_bar = 1.2 G`, and it is a four-wall channel average, the
same construction as Lau's. The superseded paragraph was the one that got
quoted.

### Confirmed from primary raw data

NASA CR-3837's appendix prints `Nu(R)`, `Nu(S)` and `Nu(AV)` per run.
`Nu(AV)` is the two-way mean of the ribbed and smooth walls to **0.084% over
33 rows**, and in a 2-ribbed / 2-smooth square duct that is exactly the
four-wall average. Reconstructing the report's printed `Gbar` column from
`St(AV)` reproduces it to +/-0.02 on 29 of 33 rows. So Han's `G_bar` is a
wall average on Han's own measurements, not on a figure label.

### What Lau's `Gbar` is actually worth

Han's `G_bar = 1.2 G` otherwise rests on a **single printed label**. Lau is
an independent lab measuring the same quantity, so committing it checks a
relationship nothing else in the dataset touches:

| | Lau `Gbar/G` | Han `G_bar/G` |
|---|---|---|
| `e+` = 106 | 1.2506 | 1.2000 |
| `e+` = 300 | 1.2410 | 1.2000 |
| `e+` = 636 | 1.2350 | 1.2000 |

Scored: **N=9, MAE 3.5%, RMSE 3.9%, bias -3.5%, 88.9% within** Lau's stated
+/-5.8% Stanton uncertainty. Han's factor sits consistently a few percent
below Lau's -- a real, recorded gap between two labs' wall-averaging ratios,
not agreement manufactured by a conversion.

### The structural tell is retired

The old argument offered a generalisable signature: Lau's ratio drifts with
`e+` (exponents 0.250 against 0.257) while Han's is constant (both 0.28), and
claimed a Prandtl normalisation must be constant while a wall-averaging ratio
need not be. The observation is correct and the inference is **falsified** --
both are wall averages, and one of them is constant. Constancy does not
identify the mechanism. Do not reuse this tell.

What does survive is the lesson item 10 already drew: a numerical agreement
at the few-percent level is not evidence of mechanism when the mechanism is
printed two pages away. This document reproduced that error rather than
learning from it, one section after citing it.

---

## Results

### `G` -- an independent confirmation

`han_1988_orthogonal` prints `G = 3.7 (e+)^0.28`; Lau prints
`4.218 (e+)^0.257`. Different coefficient *and* exponent, yet the two forms
cross near `e+ = 300` and agree across the whole overlapping range:

| `Re` | `e+` | Lau `G` | combaero `G` | |
|---|---|---|---|---|
| 10 000 | 106 | 13.98 | 13.65 | +2.4% |
| 30 000 | 318 | 18.54 | 18.57 | -0.1% |
| 60 000 | 636 | 22.16 | 22.55 | -1.7% |

Scored through the harness: **N=9, MAE 1.2%, RMSE 1.4%, bias -0.1%, 100%
within** Lau's stated +/-5.8% Stanton uncertainty.

Robust to the geometry: `han_1988_orthogonal` carries no `e/D` exponent for
`G`, and the comparison is made at matched `e+`, varying under 0.9% across
`e/D = 0.047-0.078`.

### `R` -- a ~13% disagreement that is not definitional

| `Re` | `e+` | Lau `R` | combaero `R` | |
|---|---|---|---|---|
| 10 000 | 106 | 3.618 | 3.200 | +13.1% |
| 30 000 | 318 | 3.605 | 3.200 | +12.7% |
| 60 000 | 636 | 3.597 | 3.200 | +12.4% |

Ruled out: `han_1988_orthogonal` has no `e/D` exponent for `R`, and at
`p/e = 10` its `p/e` term sits exactly at reference, so the guessed geometry
cannot be the cause. The geometry group matches for a square channel (item 9),
and both `R` are built on the ribbed-wall friction (item 5).

Converted via `R = sqrt(2/f_r) + 2.5 ln(2 e/D) + 2.5`: Han's `R = 3.2` implies
`f_r = 0.0575`, Lau's `~3.6` implies `0.0504`. **Two labs about 12% apart in
ribbed-wall friction factor for nominally the same configuration.** Recorded,
not reconciled -- see reviewer item L2.

**L2 CLOSED 2026-10-02 (#403 item 2) -- by a third lab that sides with
neither.** Taslim & Spring (1987), Northeastern / GE, measure Han's own
square-channel configuration (two opposite walls, 90 deg in-line, `P/e` 10)
and score against `han_1988_orthogonal` in `f` directly
(`taslim_spring_1987_aspect_ratio.md`). At `e/D` 0.083, the nearest they
get to Han's band, Han's friction is 14.9% LOW; on Han's basis their `R` is
2.77.

| source | lab | `R`, 90 deg, square |
|---|---|---|
| Lau (1990) | Texas A&M, Lau | ~3.6 (+12%) |
| Han (1988); CR-3837 (3.195) | Texas A&M, Han | 3.2 |
| Taslim & Spring (1987), `e/D` 0.083 | Northeastern / GE | 2.77 (-13%) |

So there is no outlier to find: the labs bracket Han by roughly +/-13% in
`R`. That spread is the honest cross-source accuracy of the friction side of
the Han family, and a fourth source would add a point to the spread, not
arbitrate it. Two caveats travel with it: Taslim's `e/D` is 6% above Han's
band, and the comparison rests on Han's own `fbar` decomposition.

Not scored by the harness: `runner.py`'s `e+` path resolves `y_axis` to `G`,
`G_bar` or `R_normalised` and has no absolute-`R` branch. The series is
committed with `scores: null` so the gap is visible rather than absent.

### Reflection folding is falsified

`han_ribbed.md` recorded reflection folding (`alpha -> 180 - alpha`) as the
bounded alternative if Eq. 4.17 were ever guarded past 90 deg. Folding
predicts V-120 behaves as V-60 and V-135 as V-45. Lau measured both:

| `e+` | `G(V-60)` | `G(V-120)` | | `G(V-45)` | `G(V-135)` | |
|---|---|---|---|---|---|---|
| 100 | 8.16 | 10.03 | **+22.9%** | 9.33 | 10.55 | **+13.1%** |
| 400 | 14.19 | 16.34 | +15.2% | 15.26 | 17.43 | +14.2% |

Higher `G` means lower heat transfer, so the reversed arrays transfer less --
and the paper's own conclusion 4 says so independently: "Reversing the 45 and
60 V-shaped rib arrays lowers the heat transfer from the channel walls and
increases the channel pressure drop."

**The open item in `han_ribbed.md` is answered: folding is not physical.** No
guard should be written on that basis.

---

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-24 | reviewer | **Third channel on Table 2, sampled as a cross.** Full first numeric column (`R` a): 3.674, 1.605, 1.232, 3.687, 1.537, 3.617, 2.070. Full last row (7/-V-135): 2.070, -0.0384, 1.992, 0.362, 3.685, 0.313. Every value matches both earlier channels exactly. Two things this buys beyond the cells themselves: the negative exponent `-0.0384` confirms the minus signs are being read -- they sit slightly detached from the digits in this scan, which was the plausible failure mode for the `R` b column -- and 0.362 lands a value in the `G` b column where OCR broke down. **Then the whole `G` b column**: 0.257, 0.355, 0.399, 0.376, 0.352, 0.324, 0.362, **and the whole `R` b column**: -0.0033, 0.0015, 0.0454, -0.0571, -0.0228, -0.0734, -0.0384. Both matching exactly. `R` b settles the sign question outright -- the column is negative, then two POSITIVE entries (0.0015, 0.0454), then four negative, and a dropped or invented minus was the plausible failure there. **Finally the remaining three columns** (`G` a, `G_bar,Lau` a, `G_bar,Lau` b), all matching. Table 2 is therefore confirmed cell by cell on three independent channels: a visual read at 300 dpi, macOS Vision OCR at 600 dpi, and the reviewer. All 42 numbers, no exceptions. That closes the one place any two channels disagreed: 0.352 and 0.324 now rest on three independent reads rather than on a visual read plus a crop, and they are the V-120 exponents the reflection-folding falsification is built from. |
| 2026-09-24 | reviewer | Channel dimensions confirmed independently: "a straight, square channel made of aluminium, 1.52 m long, cross section 7.62 by 7.62 cm". Also raised the `Gbar`/`fbar` question, which is what surfaced the false-confirmation hazard above -- it had not been checked, and the `G` comparison already made was one symbol away from being invalid. |
| 2026-09-24 | Claude | **Table 2 confirmed on a second channel: 41 of 42 numbers identical, the 42nd escalated and resolved.** `scripts/ocr_page.swift` (macOS Vision) on a 600 dpi crop reproduced the `R` a/b, `G` a, and `Gbar` a/b columns exactly, all at confidence 1.00. The `G` b column came back with only six tokens, one of them `10322` at confidence **0.30** -- the OCR merging the adjacent `0.352` and `0.324` rows and flagging itself, exactly the failure mode `extractions/README.md` describes. Escalated with a high-resolution visual crop rather than more OCR, per the same rule: the column reads `0.257, 0.355, 0.399, 0.376, 0.352, 0.324, 0.362`, confirming the original read. Those two values are the ones the reflection-folding falsification rests on, so the disagreement landed exactly where it mattered. The reviewer independently confirmed the table's structure: one case column and six numeric columns as `a, b` pairs for `R(e+) = a(e+)^b`, `G(e+, Pr) = a(e+)^b` and `G_bar,Lau(e+, Pr) = a(e+)^b`, over seven rows. That rules out a dropped row or column, and a column-pairing misread -- which is the error that would have silently swapped a coefficient for an exponent. |
| 2026-09-24 | extracted by Claude | Table 2 read off a 300 dpi render, **one channel only** -- superseded by the entry above. Definitions (Eqs. 5-8) read at 320 dpi. Geometry read from a column whose right margin the scan clips, and confirmed by the 3 inch / (3/16) inch identity. `R`/`G` findings and the reflection-folding falsification computed from the printed coefficients, not from any digitised point. |

Change **Status** at the top when reviewed, and record corrections here rather
than silently editing the tables above.
