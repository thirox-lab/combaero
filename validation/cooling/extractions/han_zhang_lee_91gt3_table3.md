# 91-GT-3 Table 3: the G_bar/G evidence

**Status: extracted, two-channel verified, 2026-09-26. Not committed as a
scored series** -- it is evidence about a CONSTANT, not data for the harness.

Source: `sources/heat_transfer/rib/han_91_GT_3.pdf`, page 9.
Han, J.C., Zhang, Y.M. and Lee, C.P. (1991), ASME 91-GT-3, "Influence of
Surface Heat Flux Ratio on Heat Transfer Augmentation in Square Channels
with Parallel, Crossed and V-Shaped Angled Ribs". Square channel,
`e/D` = 0.0625, `P/e` = 10, `L/D` = 20, `Re` 15,000-80,000.

Same three authors as figure 4.51's source (Han, Zhang & Lee 1991, JHT 113,
590-598), same rib family, so its `G_bar` is the same construction.

## Transcription (Case 1, uniform heating, is the row that matters)

`R = a(e+)^b` -- note **b = 0 for every orientation**, i.e. `R` is
independent of `e+`, which is Han's own claim elsewhere:

| | a | b |
|---|---|---|
| 90 deg rib | 3.18 | 0 |
| 60 deg crossed | 3.18 | 0 |
| 60 deg parallel | 2.05 | 0 |
| 60 deg v-shaped | 1.72 | 0 |

`G = a(e+)^b` and `Gbar = a(e+)^b`, Case 1:

| | `G` a | `G` b | `Gbar` a | `Gbar` b | `Gbar/G` at e+ 300 |
|---|---|---|---|---|---|
| 90 deg rib | 1.61 | 0.42 | 2.94 | **0.35** | 1.225 |
| 60 deg crossed | 1.82 | 0.41 | 1.90 | 0.42 | 1.105 |
| 60 deg parallel | 1.04 | 0.48 | 1.14 | 0.48 | 1.096 |
| 60 deg v-shaped | 1.01 | 0.47 | 1.63 | 0.42 | 1.213 |

Cases 2-6 (`G`) and 2-4 (`Gbar`) are in the source and were transcribed but
are not reproduced here: they vary the ribbed-to-smooth wall heat flux
ratio, which no combaero element expresses.

## Verification

Visual read at 400 dpi plus macOS Vision OCR at 400 dpi. The `R` and `G`
sub-tables agree cell for cell on both channels (one OCR `0.S1` for 0.51,
unambiguous).

**One gap: OCR never captured the `Gbar` 90 deg `b` column at all.** That
column also holds **0.35**, the only value in the entire `Gbar` table
outside 0.39-0.48 -- unverified and an outlier at once. Escalated to a
900 dpi crop per `extractions/README.md` (a crop, not more OCR) and it
reads cleanly: 0.35 / 0.42 / 0.40 / 0.42. Reviewer asked to check it.

## What it was extracted for, and what it settled

To answer whether `G_bar/G` generalises well enough to replace Han's
published constant 1.2 with a correlation. **It does not.**

Pooling every closed-form measurement available -- this table, Lau's table
2, and CR-3837's per-run `Nu(R)`/`Nu(AV)` split -- gives 16 configurations
across three rigs at `e+` = 300:

    min 1.096   max 1.413   mean 1.284   spread 29%

and within a configuration it drifts with `e+` too (this table's 90 deg rib
runs 1.29 down to 1.13 over `e+` 150-1000, so even at 90 degrees Han's
printed 1.2 is a mid-range approximation of a moving quantity).

**69% of that variance is BETWEEN RIGS, not within them:**

| source | n | mean | internal spread |
|---|---|---|---|
| 91-GT-3 | 4 | 1.160 | 12% |
| CR-3837 | 5 | 1.323 | 3% |
| lau1990 | 7 | 1.327 | 14% |

### "Three rigs" is not three laboratories

91-GT-3 and CR-3837 are Han's own laboratory; `lau1990` is Lau's, in the
same department, and Lau sat on the committees for the others. The phrase
reads as more independent than it is, and the consequence is two-sided:

* **the "no generalised correlation" conclusion is STRONGER.** Three rigs
  inside one laboratory tradition differing by 69% of the variance
  between them makes the offset rig-level -- geometry, instrumentation,
  data reduction -- rather than an artefact of comparing institutions.
  There is even less for an angle or shape term to reach.
* **the "MAE 8.3%" is OPTIMISTIC.** This is not a sample of the world's
  rigs. The 29% spread is a LOWER BOUND on the true one, and Han's 1.2
  may score worse against a non-TAMU rig than it does against this pool.

Same monoculture as #403 item 2, where it blocks the Lau friction
disagreement outright. Noted 2026-10-01.

No angle or shape term can reach a rig offset, so a "generalised" `G_bar/G`
would be fitting rig identity. Against the pooled population Han's 1.2
scores bias -6.1%, MAE 8.3%.

**Decision: keep 1.2, apply it as published at every rib angle, and report
what it costs as accuracy.** The runner previously refused off-90 series,
which withheld a number Han does publish and turned a known accuracy limit
into a missing answer. Pinned by
`test_keeping_hans_constant_is_justified_by_the_pooled_measurements`, which
fails if the between-rig share drops -- the condition under which a
generalised correlation would become worth revisiting.

## Not used

- Cases 2-6: wall heat flux ratio has no expression in combaero.
- The `R` coefficients: 60 deg parallel 2.05 and v-shaped 1.72 are a
  SHAPE-dependent friction result that `han_park_1988_angled` cannot carry
  (it has angle only). Relevant to reviewer item L1, not to this decision.
