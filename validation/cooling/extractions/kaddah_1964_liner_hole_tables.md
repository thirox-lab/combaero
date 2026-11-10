# Kaddah (1964): liner air-entry hole Cd and jet angle under supply crossflow

**Source.** Kaddah, K. Sh. (1964). *Discharge coefficient and jet deflection
studies for combustor liner air-entry holes.* Diploma of Advanced Engineering
thesis, College of Aeronautics, Cranfield. Vol. 1, Tables 11-30 (PDF pages
99-118). Vol. 2, Figs. 35-103, re-plots them. Both volumes are registered in
`validation/sources.yaml`; the scans stay in `sources/`.

**Data.** `validation/cooling/data/kaddah1964/`:

| file | what |
|---|---|
| `table11_C1.csv` ... `table30_R8.csv` | the 20 printed tables, **as printed**, in the thesis's own columns and units (lb/hr, lb/cu.ft, ft/sec, deg) |
| `tables.yaml` | rig, definitions, per-table geometry, a*, t/d, the Vol. 2 figure numbers, the liner-envelope class, and the table-level notes |
| `cell_flags.csv` | per-cell read confidence and printed misprints, each with its alternative |

The tables are multi-column printed tables, not (x, y) figure series. So they
do not use `schema.py`'s `metadata.yaml` (and `load_dataset()` does not read
them); the eriksen1971 as-printed tables are the precedent.
`kaddah_runner.py` reads them. The only change from the page: a decimal comma
or semicolon the typist used is written as a dot.

## Rig and definitions (Vol. 1, Sections 3, 7 and Appendix 2)

- **Duct and plate.** 8 in wide x 1.25 in deep, 18 in long, so A1 = 10 sq.in.
  1/16 in aluminium plate, square-edged holes flush with the duct floor.
- **Pressures.**
  - p1: approach static at station 2, 7.875 in upstream of the hole centre, held at 33 in Hg abs.
  - pj: the plenum static.
  - P1: pitot total at station 3.
- **Configurations.** Circular C1-C5 (d 0.329-1.59 in, t/d 0.19-0.039). Oval O1-O7 and rectangular R1-R8,
  aspect ratio 1-8. a* = n Ah/A1 runs 0.0595-0.604.
- **C.** `C = Lj/Lth`, with:
  - `Lth = 25 rhoj Vj n Ah` (lb/hr, Ah in sq.in);
  - `rhoj = rho1 (pj/p1)^(1/1.4)`;
  - `Vj = 18.28 sqrt(h/rhoj)` (A2.1.1-A2.1.4).

  The printed Vj is reproduced by `sqrt(2 g (P1 - pj)/rhoj)` at p1 = 33 in Hg (median +0.2%).
  So h = P1 - pj, and **C is referenced to the approach TOTAL pressure** (`conventions.py`: `Cd_total`).
- **K.** `K = (P1 - pj)/(P1 - p1)`, computed as `(rhoj/rho1) Vj^2/V1^2`. This is Rohde's
  velocity head ratio and TN 3663's Phi.
- **Other columns.**
  - `alpha = Lj/L1` (percent), and `mu = alpha/a*`.
  - phi is the angle between the jet axis and the downstream plane of the hole; 90 deg is a normal jet,
    the same sense as combaero's `injection_angle_deg`.
  - delta needs the hole-centre static, which is not tabulated.
- **No-cross-flow blocks.** The duct is closed downstream (alpha = 100%), and V1 is the residual feed.

## Transcription (user sign-off 2026-11-10)

1. **Two independent channels.**
   - macOS Vision OCR with word boxes, grouped into rows.
   - A visual read of the page images by two readers who never saw the OCR.

   They agreed on 3325 of 3766 cells. The 441 disagreements resolved as follows:
   - 303 were OCR noise (the decimal point read as a digit).
   - 78 were settled by the table's own arithmetic or by physical plausibility.
   - 53 were settled by a third read of the scan.
   - 7 were left open and taken as the third reader's pick (the user's decision). They are flagged `sequence_judged`.
2. **The thesis's own arithmetic, every row.** The checks are C = Lj/Lth, alpha = Lj/L1, mu = alpha/a*,
   Lth = 25 rhoj Vj nAh, the isentropic rhoj, V1 = L1/(rho1 A1), K, and Vj from P1 - pj.
   The median of every check is within +/-1% (K: -1%).
3. **Cross-check against Vol. 2.** Each table was overlaid on its C and phi vs p1/pj figure (Figs. 35-39,
   49-53, 56-57, 71-75, 78-80); the axes were calibrated from OCR'd tick labels. The overlays agree within
   plotting precision, with two exceptions:
   - **Fig. 35 (C1)** plots the no-cross-flow C 0.03-0.07 below Table 11, and two 130 ft/s points
     differently. The table is primary: it is self-consistent, C = Lj/Lth to 0.1%.
   - **Fig. 49 (O1)** supports .602 for the misprinted Table 16 130/4 C.

`cell_flags.csv` holds the per-cell read confidence:

| flag | n | meaning | scoring |
|---|---|---|---|
| `sequence_judged` | 43 | ambiguous digit (overstruck, faint, overtyped, handwritten) judged from its shape and the row sequence; NOT exact | rows whose Cd inputs carry it are reported separately |
| `arithmetic_confirmed` | 36 | ambiguous digit; the table's arithmetic holds for the chosen read and fails for every alternative | normal |
| `printed_inconsistent` | 7 | the printed value contradicts the table's own arithmetic; kept as printed, alternative recorded | excluded by default |

These are the `printed_inconsistent` cells:

| table | cell | printed | alternative |
|---|---|---|---|
| 11 | no-cross-flow row 4, phi | 840 | 84 |
| 13 | 80 ft/s row 1, Lj | 179 (last digit illegible) | 197 |
| 14 | 80 ft/s row 2, C | .488 | .448 |
| 16 | 130 ft/s row 4, C | .502 | .602 |
| 17 | 130 ft/s row 1, C | .435 | .454 |
| 24 | 80 ft/s row 3, alpha | 18.6 | 28.6 |
| 28 | 130 ft/s row 3, Lth | 1590 | 1960 |

Two problems are table-level, not cell-level:
- **Table 27 (R5).** Every mu and Lth implies n Ah = 6.72 sq.in (a* 0.672), not the 6.0 sq.in (a* 0.600) of
  Section 7.3 and Fig. 75. If the slots are 6.0 sq.in, the printed C are about 11% low. The table is
  excluded by default.
- **Table 13 (C3).** mu implies a* 0.246 against the stated 0.2542. mu is not a scoring input.

Table 15 (C5) is handwritten.

## Conventions used for scoring

- **Cd basis.** Static, via the registered `Cd_total -> Cd_static`, `sqrt(K/(K-1))`. The relative error is
  identical on the total basis, because one factor multiplies both sides.
- **U1/Vi, both definitions, always labelled.** Incompressibly, Vi/V1 = sqrt(K - 1), with Vi from the duct static.
  - Kaddah's **approach** V1: `U1/Vi = 1/sqrt(K - 1)`.
  - combaero's **mean** of the approaching and departing flow (`crossflow_velocity_ratio`):
    `U1/Vi = (1 - alpha/2)/sqrt(K - 1)`.

  Kaddah's bleed runs from 2.7% to 100%, so the two definitions differ by up to 2x. In Rohde and in TN 3663
  they coincided.
- **Liner envelope (score first).**
  - **Inside:** C1, O1, R1 (a* 0.0595) and C2, O2, R2 (a* about 0.12).
  - **Marginal:** C3, O3, R3, R6 and O6 (a* about 0.25, slot length no longer than the duct depth).
  - **Outside:** the rest (a* 0.38-0.60, or slots longer than 1.25 in).
- **Confound.** C1-C5 vary d at a fixed t, so t/d and a* change together; a wall-thickness term cannot be
  isolated in this set.

## Scores -- ACCURACY, CROSS-SOURCE

`uv run python -m validation.cooling.kaddah_runner` prints the full report: modes, U1 definitions,
envelope classes, read confidence, and K bands. Each score is predicted/measured - 1 on Cd_static.

**Headline: U1 = mean** (combaero's mean of the approaching and departing duct flow, which is what a user
of the element gets); the approach-velocity rows (Kaddah's own V1) are kept beside it and are never
dropped (user decision, 2026-11-10). The runner itself takes no default.

**Eq. 17** is anchored on each table's own no-cross-flow block (the source's own way of validating Eq. 17)
and scored on the 80 and 130 ft/s rows. Clean inputs:

| U1 | envelope | n | bias | MAE | K 1-2 | K 2-10 | K >= 10 |
|---|---|---|---|---|---|---|---|
| approach | inside | 54 | -5.3% | 7.1% | -32.8% (3) | -5.3% (25) | -2.2% (26) |
| mean | inside | 54 | -3.0% | 5.8% | -30.1% (3) | -2.6% (25) | -0.2% (26) |
| approach | marginal | 36 | -13.6% | 13.6% | -40.2% (4) | -14.4% (18) | -5.0% (14) |
| mean | marginal | 36 | -12.0% | 12.0% | -39.8% (4) | -12.2% (18) | -3.7% (14) |
| approach | outside | 68 | -16.0% | 16.5% | -35.7% (14) | -14.2% (37) | -3.9% (17) |
| mean | outside | 68 | -15.8% | 16.3% | -38.6% (14) | -13.1% (37) | -2.9% (17) |

The rows with sequence-judged inputs are reported separately: 2 inside, 5 marginal and 2 outside.

Inside, Eq. 17 under the approach definition, by configuration (80 / 130 ft/s bias):

| config | 80 ft/s | 130 ft/s |
|---|---|---|
| C1 | +6.6% | +4.2% |
| O1 | -8.0% | -12.1% |
| R1 | -1.0% | -14.3% |
| C2 | -4.2% | -5.4% |
| O2 | -2.8% | -7.1% |
| R2 | -7.2% | -10.3% |

C1 is anchored on the no-cross-flow block that Fig. 35 disputes; a lower anchor would lower C1's predictions.

**The full chain** (`mcgreehan_schotsch_1988_cd`, r/d 0, L/d = t/d, Re from rhoj Vj d/mu) runs on circular
holes only (the slots are refused) and scores every row, the no-cross-flow baseline included. Clean inputs:

| U1 | envelope | n | bias | MAE |
|---|---|---|---|---|
| approach | inside (C1, C2) | 26 | -7.3% | 7.9% |
| mean | inside (C1, C2) | 26 | -6.7% | 7.6% |
| approach | marginal (C3) | 13 | -11.0% | 11.0% |
| mean | marginal (C3) | 13 | -8.6% | 8.6% |

The chain's no-cross-flow baseline reads 3-7% low on C1-C3.

**Verdict (#807).** At liner geometry and moderate crossflow (K >= 2), Eq. 17 is within a few percent. Under
strong crossflow (K < 2) it reads **30-45% LOW** in every class and under both U1 definitions. Thin liner
holes keep more Cd than Eq. 17 predicts. TN 3663 (NACA Lewis) showed the same direction, so this is now
seen by two independent labs. Rohde's thick hole reads the opposite way. The low-K rows are few at small
a*: 3 inside, 4 marginal, 14 outside.

## Jet angle -- REPORTED, NOT SCORED (#813: no closure yet)

These are the 80 and 130 ft/s rows, error = predicted - measured in degrees:

| envelope | n | Report 1373 Eq. 9 (cos phi = V1/Vj) | NREC Eq. 97, Cd_inf 0.625 | NREC Eq. 97, Cd_inf = own no-cross-flow C |
|---|---|---|---|---|
| inside | 57 | +0.7 / 3.4 | +8.0 / 8.8 | +4.5 / 6.3 |
| marginal | 42 | -2.6 / 2.8 | +5.7 / 7.5 | +6.3 / 7.1 |
| outside | 71 | -5.8 / 6.3 | +2.3 / 5.7 | +5.7 / 8.4 |

Each entry is bias / MAE. At liner scale, Eq. 9 ("the jets keep the approach velocity") is the better of
the two printed forms.
