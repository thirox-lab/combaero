# NASA CR-3837 (Han, Park and Lei, 1984) -- tabulated primary

**Status: extracted and committed, 2026-09-26.** Two-channel verified.
Source: `sources/heat_transfer/rib/NASA_CR-3837_1984.pdf` (gitignored), 164 pp,
grant NAG3-311, November 1984. Journal version Han, Park and Lei (1985),
ASME J. Engng Gas Turbines Power **107**, 628-635.

## Why this source matters

It is the first source in the cooling dataset whose **every run is
tabulated**. Appendix 7.3 prints, per run:

```
Re | m | Tw(R) | Tw(S) | Tb | q"(R) | q"(S) | Nu(R) | Nu(S) | Nu(AV) | St(AV) | e+ | St/St | f/f | n | P/P | Rbar | Gbar
```

The `Nu(R)` / `Nu(S)` split is unique here. Everywhere else the ribbed-wall
Stanton number has to be assumed; here it is read.

## Scope

Square duct only (`W/H` = 1), `e/D` = 0.063 only, `P/e` = 10 and 20, rib
angles 90/75/60/45/30/15, two entrance conditions. Committed: the
long-entrance heat-transfer runs at `P/e` = 10, page 128, `X/D` = 11.5 --
33 runs across five angles.

## The circularity gate (#391)

| against | verdict |
|---|---|
| `han_park_1988_angled` (Eq. 4.17/4.18) | **independent** -- scored |
| `han_1988_orthogonal` (R = 3.2, G = 3.7(e+)^0.28) | **circular** -- never score |

Han and Park (1988) measured `e/D` = 0.047 and 0.078 only, on contract
NAS3-24227, with its own data report CR-4015 (its ref [9]). Its conclusion
(4) says the square-channel results "confirm observations in the previous
study [8]" -- CR-3837's journal version -- i.e. compared against, not
pooled into. Eq. 4.17/4.18's validity range 0.047-0.078 brackets 0.063
without having measured it.

Figure 4.46's right box is a different story: it prints "Han (1984),
8,000 <= Re <= 80,000" at `e/D` 0.063 for both `P/e`, a four-way fingerprint
match to CR-3837's Table 2. Those marks are these runs, and figure 4.46 is
where `R = 3.2` and `G = 3.7(e+)^0.28` were fitted.

## The conversion, and why it is not optional

The appendix prints `Rbar` on the **channel-average** friction `fbar`, and
`Gbar` on the **four-wall average** Stanton number. combaero uses Han's
four-sided equivalent `f` and the **ribbed-wall** Stanton number.

The printed `Rbar` is **5.35** where the same measurement on Han's basis is
**3.195** -- a 67% gap from convention alone. Committing the printed columns
would have been a silent definitional error of exactly the kind #389 names.

```
fbar = 2 / (Rbar - 2.5 ln(2 e/D) - 2.5)^2
f_s  = 0.046 Re^-0.2          reproduces CR-3837's tabulated smooth duct to ~1%
f    = 2 fbar - f_s           Han's four-sided equivalent at W/H = 1
St_r = St(AV) Nu(R)/Nu(AV)    ribbed wall, from the printed split
e+   = (e/D) Re sqrt(f/2)
R    = sqrt(2/f) + 2.5 ln(2 e/D) + 2.5
G    = R + (f/(2 St_r) - 1) / sqrt(f/2)
```

## Two-channel verification

Page 128, all five angle blocks, on a visual read at 200-400 dpi and macOS
Vision OCR at 400 dpi upscaled. **Every discrepancy was an OCR error caught
by an independent check**, and no cell survives where the channels disagree:

| block | OCR said | settled by |
|---|---|---|
| 90 | `Nu(AV)` 247.6 | `(62.6 + 32.6)/2 = 47.6`, the mean identity |
| 75 | `Gbar` 22.57 | reconstruction gives 21.59 |
| 75 | a garbled duplicate row | OCR noise, no counterpart on the page |
| 30 | `Re` 152,257 | outside the stated 7,000-90,000; border stroke read as "1" |
| 30 | `Rbar` 1.58 | breaks the monotone 8.17 -> 5.93 sequence |

OCR also **recovered a 90 deg row the visual pass had dropped** (Re 16,544),
which prints no `e+` and was skipped as incomplete. It is committed.

Blocks 60 and 45 agreed with no exceptions; OCR dropped a few cells but
contradicted none.

## Result

Out-of-sample against Eq. 4.17/4.18, scored through the harness:

| alpha | `G` MAE | `G` bias | within 6.8% |
|---|---|---|---|
| 90 | 8.4% | -8.1% | 42.9% |
| 75 | 3.7% | -0.8% | 83.3% |
| 60 | 11.6% | -11.6% | 16.7% |
| 45 | 11.6% | -10.1% | 14.3% |
| 30 | 19.9% | -19.9% | 0.0% |

**`R` (Eq. 4.17) is confirmed, and now scored through the harness** since
the absolute-`R` branch landed (#339 plan 1a):

| alpha | `R` MAE | `R` bias | within 6.6% |
|---|---|---|---|
| 90 | 3.9% | -2.3% | **100%** |
| 75 | 7.6% | -3.1% | 83.3% |
| 60 | 6.1% | +6.1% | 83.3% |
| 45 | 1.9% | +1.7% | 85.7% |
| 30 | 13.7% | +10.4% | 57.1% |

Four of five angles sit inside the source's own 6.6% friction uncertainty
for most of their points. Only 30 deg is poor, which is also where Eq. 4.18
reads furthest low -- worth noting that the shallowest angle is the weakest
for both quantities, and not explained here.

**`G` (Eq. 4.18) reads systematically low**, and the bias is negative at
every angle. That bears on **item 23** of `han_ribbed.md`, an unresolved
internal disagreement between figure 4.46 and Eq. 4.18 of up to 22% with
different `e+` exponents. An independent dataset now lands on the figure's
side. Recorded, not reconciled.

## What this extraction found in the harness

Scoring these series surfaced a defect in `runner.py` that predates them:
on the `e+` path, `run_series` copied `e_D`, `p_e` and `W_H` out of
`series.geometry` but **never `alpha_deg`**, which is a top-level field.
Every angled series was scored at the probe default of 90 degrees --
against `han_park_1988_angled`, the one set whose entire subject is rib
angle. Fixed, with a test. See the review log in `han_ribbed.md`.

## Stated uncertainty

Average friction factor < 6.6%, Nusselt number < 6.8%, both "for Reynolds
number greater than 10,000". Five committed runs sit below that (Re 8,314 /
6,281 / 6,385 / 6,571 / 9,637) and the source makes no claim for them.

`St_r` is on **projected** area -- "not including the increased rib surface
area", the source's own words, matching Lau's convention.

## Not done

- `P/e` = 20 runs, and the sudden-contraction entrance block, are not
  extracted. Both are in the same appendix.
- `alpha` = 15 is not committed: Eq. 4.17/4.18's validity stops at 30.
- The friction appendix (7.2) is not extracted; `f_s` comes from the
  Blasius form Han and Park themselves use, checked against the tabulated
  smooth duct to ~1% rather than read run by run.
