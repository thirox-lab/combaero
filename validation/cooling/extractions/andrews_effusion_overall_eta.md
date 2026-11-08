# Effusion overall cooling effectiveness -- Andrews 88-GT-290 Figure 10

Andrews, G.E., Asere, A.A., Hussain, C.I., Mkpadi, M.C. and Nazari, A.
(1988). "Impingement/Effusion Cooling: Overall Wall Heat Transfer."
ASME 88-GT-290, Gas Turbine and Aeroengine Congress, Amsterdam.
`sources/heat_transfer/film/Andrews_88-GT-290.pdf`

This closes #387's last piece -- the overall effectiveness OUTPUT -- and
reshapes it. The internal side is `andrews_effusion_internal_h.md`; the
external side is `baldauf_2002_film_effectiveness.md` and
`film_superposition.md`.

Later addition (2026-11-08): Behrendt et al. (2008), GT2008-51191, is a
reacting 5-10 bar heated-wall `eta_tot` set (effusion, 20 and 30 deg holes;
coolant temperature, pressure drop and pressure varied). It is data only,
unscored, and a cross-source accuracy set for the coupled wall, not a
gas-side closure. See
[behrendt_2008_effusion_eta_tot.md](behrendt_2008_effusion_eta_tot.md).

## What the figure is

Page 7, `eta = (T_g - T_w)/(T_g - T_c)` (the paper's Eq. 3) against
coolant mass flow per unit plate area. Eight curves: effusion B and C
alone, impingement A alone, the A/B and A/C combinations, and three
other technologies for comparison (porous wall, Lamilloy, Transply).

Only B and C are scored. A/B, A/C and A alone need an impingement closure
that does not exist here; the other three are other people's measurements
of other technologies.

## The rig, as the paper describes it

> "the wall of a 76 mm by 152 mm wide air cooled duct through which the
> product gases from a propane preheater flowed at 750K and a Mach number
> of 0.05" -- page 3

That gives `D_h = 101.3 mm`, `U_g = 26.8 m/s`, `Re_Dh = 38375` and a
Dittus-Boelter `h_0 = 46.1 W/m2K`. The paper's own statement of a
"coolant to hot gas density ratio of approximately 2.5" falls out of the
two temperatures at 2.54, which is a consistency check on the conditions
rather than an extra input.

Geometry is Table 1: plate B `D = 2.16 mm, X/D = 7.1, A/A_h = 5.28`,
plate C `D = 3.27, X/D = 4.7, A/A_h = 3.4`, both 6.3 mm Nimonic 75. The
printed `A/A_h` values reproduce from `(X^2 - pi D^2/4)/(pi D L)` to 1%,
so the table is self-consistent with the stated thickness.

## THE FIGURE AND THE TEXT DISAGREE, TWICE

Recorded, not resolved:

| | the plot prints | the body text says |
|---|---|---|
| temperatures | `Tg = 744K, Tc = 293K` | "a Tg of 750K and Tc of 295K" (p. 7) |
| gap | `Z = 6.4mm` | "An 8 mm impingement gap was used throughout" (p. 3) |

`Z` is defined as the impingement gap, so the second is a real conflict --
but it does not touch the effusion-only curves, which are the scored ones.
The first is worth under 0.5% on `eta`, which is a ratio; the text's
values are used and the difference is reported.

## The closure, and what it lumps

    eta = h_i / (h_i + h_gas)

`h_gas` is supplied by the caller. There is no separate film term, and
the reason is NOT that there is no film.

### An adiabatic effectiveness is half of a pair

The two-temperature form an external film really needs is

    q = h_f (T_aw - T_w),   T_aw = T_g - eta_f (T_g - T_c)
    eta = (h_i + h_f eta_f) / (h_i + h_f)

and it takes both `eta_f` and `h_f`. Baldauf supplies only `eta_f`, and
structurally cannot supply `h_f`: an adiabatic wall passes no heat, so
the experiment fixes `T_aw` and says nothing about the coefficient.

### What Baldauf DOES contain -- an earlier reading of this had it backwards

The counter-rotating vortex pair entraining hot gas under a lifted jet
changes `T_aw`, which is exactly what IR thermography on an insulated
wall records. **Baldauf has that**, in Eq. (38)'s decay branch, and
steeply. The falling exponent is a function of hole spacing alone:

| | Eq. | plate B (s/D 7.06) | plate C (s/D 4.66) |
|---|---|---|---|
| decay exponent `b_pk` | (12) | 4.57 | 3.24 |
| peak location `mu_0` | (14) | 2.25 | 1.13 |
| peak height `eta_c0` | (15) | 0.137 | 0.228 |

Measured off the shipped correlation at `x/D = 10`, the net large-`M`
slope is -3.0 (B) and -3.2 (C), plate C matching Eq. (12) almost exactly.

Andrews' own sentence -- "the extent of the jet stirring of the cooling
film determined the effectiveness of the film cooling. This stirring of
the boundary layer was reduced as the hole size was increased due to the
lower jet velocities at a fixed mass flow" -- describes film DESTRUCTION.
That is the mechanism Baldauf models, not a missing one. An earlier draft
of this record attributed it to the missing coefficient; it is corrected
here because the two are different physics and conflating them points the
next reader at the wrong fix.

### The measurement, read both ways

    eta_f admissible at h_f = h_0 = (eta_meas (h_i + h_0) - h_i)/h_0

| | offered alone at the smooth-duct h_0 | admitted with its augmentation |
|---|---|---|
| plate B | 0 of 66 points admit any; max -0.170 | needs `F = h_f/h_0` 2.7 to 4.7 |
| plate C | 10 of 72 admit some, all `G < 0.24`; max +0.090 | needs `F` 1.6 to 3.3 |

against 0.27 to 0.58 from Baldauf superposed over the plate's ten rows.
Only the PAIR is identifiable from one eta curve, so the caller supplies
it as one `h_gas`. That form has zero free parameters, which is what
makes it falsifiable where the split form fits anything.

**Plate C's +3.1% is a cancellation, not an absence.** A real film raises
`eta` and its augmentation lowers it, and at that geometry the two nearly
balance. Saying so matters: the agreement is not evidence that the
closure is right in general.

### The absolute level is confounded; the ratio is not

The test plate sits only `x/D_h = 1.5` into the duct, deep in the entry
region, where Dittus-Boelter's fully developed value understates `h_0`. A
standard entry correction raises it 1.76x -- 46.6 to 82.0 W/m2K -- and
brings plate C's required `F` to 0.9-1.9, ordinary film-cooling augmentation.

### The augmentation each plate requires, against each baseline (#476)

`F = h_f / h_0` from `implied_gas_side_h`, every point; recomputed 2026-10-20 on
main dfae3de; the ranges quoted before were stale. This table owns these numbers; the element
docstring, `docs/API_PYTHON.md` and `conventions.py` point here.

| plate | film admitted | n | vs Dittus-Boelter (46.6 W/m2K) | vs entry-corrected (82.0 W/m2K) |
|---|---|---|---|---|
| B | none | 66 | 1.33-2.71 | 0.76-1.54 |
| B | Baldauf, Sellers over 10 rows | 66 | 2.69-4.68 | 1.53-2.66 |
| C | none | 72 | 0.82-1.36 | 0.47-0.77 |
| C | Baldauf, Sellers over 10 rows | 72 | 1.57-3.30 | 0.89-1.88 |

A third baseline, the turbulent flat-plate boundary layer at the panel's
streamwise position that published `h_f/h_0` ratios use, is now offered by the
element (`gas_baseline="flat_plate"`, CR-2648's St_o, with `gas_start_length`
the panel's x from where the boundary layer starts). It is NOT applied to
this rig: 88-GT-290 does not print where the hot-side boundary layer starts
relative to the plate (only that eta was taken 127 mm downstream of the
impingement gap's leading edge, p. 3), so any start length would be invented.
It is scored on Behrendt 2008, whose bleed fixes it
(extractions/behrendt_2008_effusion_eta_tot.md).

Not applied: the runner reports the rig as the paper describes it.
Recorded so the absolute `F` values are not over-read. It divides out of
every ratio below.

## Result

| | bias | MAE | n |
|---|---|---|---|
| plate C (X/D 4.7) | **+4.0%** | 5.0% | 72 |
| plate B (X/D 7.1) | **+24.6%** | 24.6% | 66 |
| basis | **accuracy** (1986 correlations, 1988 data) | | |

(re-scored 2026-10-08 after #485: air's viscosity was 11.8% low before); it was +3.1% / +23.7%.

Reported as two rows, never pooled: a pooled +13% would describe neither.
The runner's reporting group is per plate so the rollup cannot form it.

**Per plate and NOT per jet regime.** Splitting the points by velocity
ratio looks informative -- roughly +7% below `VR = 1` against +19% above
-- but both plates cross `VR = 1` inside their own `G` range, so the
split reports the plate mix:

| | plate B | plate C |
|---|---|---|
| VR < 1 | +21.2% (n=17) | +2.2% (n=52) |
| VR >= 1 | +24.6% (n=49) | +5.6% (n=20) |

The plate moves the error by a factor of ten; the regime label by a few
points.

## THE FINDING: a gas-side residual no film correlation can close

Plate B needs about 1.7x plate C's gas-side coefficient, and `h_0`
divides out of that so it rests on no assumption about the rig:

    B/C required h_gas      G=0.6   G=1.0   G=1.4
      no film term           1.82    1.85    1.96
      Baldauf's own film     1.75    1.66    1.79

Scaling Baldauf's `eta_f` from 0 to 1.25x moves the absolute `F` a long
way -- plate B's at `G = 0.6` runs 2.07 to 5.99 -- and the ratio from
1.96 to 1.60. **So the split is bounded away from 1 whatever film model
is chosen.** A better film correlation cannot close it, because the
missing quantity is a coefficient and an adiabatic measurement cannot
yield one.

Pinned by `test_no_film_correlation_of_any_magnitude_closes_the_split`.

### No threshold is offered

The enhancement was tested for collapse against the velocity ratio, the
blowing ratio and the momentum flux ratio. It collapses on none: at equal
`VR` the plates still differ by 1.5-1.6x, at equal `M` by 1.35-1.45x, at
equal `I` by 1.3-1.5x. Two plates cannot separate a jet parameter from a
geometry one, so `jet_regime()` reports all three and fits nothing.

Jet velocities, for orientation: plate B ejects at 53 m/s into a 26.8 m/s
crossflow at `G = 1.0`; plate C at 23 m/s, below it.

## Why no coolant heat-up term

The coolant gains 12 K at `G = 1.4` and 60 K at `G = 0.2` passing through
the wall, so an explicit heat-up looks obviously missing. It is not:
86-GT-225's `h` is fitted from the plate's transient cooling rate against
the SUPPLY temperature (Eqs. 1-2 of 88-GT-290), so the heat-up is already
inside it.

The evidence is that adding it does not FLATTEN the error against `G` --
a genuinely missing `G`-dependent term would. It tilts it: plate C runs
-2.8% to +5.8% without (8.6 points of span) and -8.8% to +4.2% with
(13 points). Wall conduction is omitted for the same reason -- the
transient technique lumps the plate -- and is 3.5% of `1/h_i` anyway.

## Digitisation

Eight curves, digitised by the user from page 7.

**The page is rotated and the raw picks are committed anyway.** The frame
corners give 0.94% skew read as y-per-x-span and 0.96% read as
x-per-y-span -- two independent measurements of one real rotation. But
the digitiser's own affine already absorbed it: applying the four-corner
inverse makes the y ticks WORSE, from a worst error of 0.00178 (0.30% of
span) to 0.00255, and leaves the x ticks unchanged at 0.44%. The third
channel below cannot separate them either (0.57 against 0.52 percentage
points). Uncorrected, on the evidence.

This is the Rohde Fig. 10 check run again and coming out the other way,
which is the reason to run it rather than assume.

**Table 3 is a third channel and it validates the curve LABELS.** Page 7
tabulates the ratios BETWEEN these curves at `G = 0.2, 0.5, 1.0` under
"Relative Cooling Effectiveness", without saying what the ratio is.
Reading it as `(eta_num - eta_den)/eta_den` reproduces all nine evaluable
printed cells to 0.57 percentage points mean, worst 2.0:

| G | A/B vs A_imp | A/C vs A_imp | A/B vs B | A/C vs C |
|---|---|---|---|---|
| 0.5 | 6.5 vs 6 | 17.0 vs 17 | 26.0 vs 24 | 20.0 vs 20 |
| 1.0 | 5.2 vs 5 | 12.2 vs 13 | 26.0 vs 26 | 18.1 vs 19 |

Nine cells agreeing is not chance, so the READING is confirmed as well as
the digitisation. The table's dashes in the `A_imp` columns at `G = 0.2`
match the A_imp curve starting at `G = 0.487` -- an independent
confirmation that the curve identification is right. One printed cell
(A/B vs B at `G = 0.2`) is not evaluable: the A/B curve's first point is
at `G = 0.203`, and inventing a point to fill a validation cell is not
done.

## Candidate source for the gap -- NOT obtained, NOT verified

**Badal et al.**, a three-regime film-cooling correlation: fully attached
(no jet lift-off), fully lifted, and an interpolated transition, selected
on blowing ratio. IR rig on a scaled nozzle guide vane geometry,
cylindrical and fan-shaped holes. Known to this record only from an
abstract supplied by the user; the paper has not been read and nothing
here rests on it.

It would address the `eta_f` half in the lift-off regime, which is a real
gap -- Baldauf is extrapolating hard here (`M` to 7.0 against a 2.5
limit, density ratio 2.54 against 1.8, plate B's `s/D` 7.06 against 5).

**But it would not close THIS gap.** The bounding result above says no
choice of `eta_f` brings the two plates together, and an IR adiabatic
measurement yields `eta_f` and not `h_f`. What this needs is a gas-side
augmentation `h_f/h_0` for full-coverage effusion -- a different
measurement, on a heated wall.

## Deliberately not implemented

**A fitted gas-side augmentation.** Doubling `h_g` takes plate B from
+24.6% to +2.7%, so a runner that "calibrated" it would report both
plates as good fits and erase the finding entirely. That is why
`gas_side_sensitivity()` reports the multiplier and never applies it.

**An overall-effectiveness correlation as the closure.** It would already
contain the internal convection the network computes, and the two would
double-count. `eta` is an output here and must stay one.

## The knobs, and that they reach the target

The library does not correct the 25% miss on plate B -- matching a rig is
the user's job. But a knob that is offered and does not reach is worse
than no knob, so each was verified against plate B at `G = 0.6`
(measured 0.6074, untuned 0.7665, +26.2%; values after #485):

| knob | value needed | result |
|---|---|---|
| `gas_augmentation` | 2.123 | exact to 1e-16 |
| `eta_film` 0.318 + `gas_augmentation` 4.454 | the physical pair | exact |
| `internal_Nu_multiplier` | **0.471** | exact -- and WRONG |

**Both single knobs reach the target; only one is attributable.** The
internal multiplier gets there by halving the coolant-side coefficient,
which is the wrong direction on the evidence: against Andrews' own Fig. 8
that correlation runs 7.1% LOW, so correcting it would RAISE `h_i`. The
number it needs is itself the evidence that it is the wrong dial. Using
it to absorb a gas-side error would be reward hacking with a user-facing
knob, and `test_the_internal_knob_also_reaches_it_but_should_not_be_used`
records that.

**One scalar is a rig match, not a model.** Fitting `gas_augmentation` at
`G = 1.0` takes plate B from 24.6% MAE to 5.7% -- a real improvement --
but the residual runs **-28.4% to +3.5%** and is strongly asymmetric,
because the required augmentation FALLS at low G (1.3 at `G = 0.15`
against 2.4 at 1.0). The average closing hides the low-G tail.

`internal_Nu_multiplier` was added by this work; the element previously
had no coolant-side tuner at all, which the verification exposed.

## Falsification

Eight perturbations, each applied alone and reverted.

| perturbation | result |
|---|---|
| `predict()` adds Baldauf's film back in | RED x3 (both plates, element/runner tie, rollup rows) |
| reporting group per jet regime instead of per plate | RED x2 |
| gas side uses the duct WIDTH, not the hydraulic diameter | RED |
| the coolant heat-up term is added to `predict()` | RED x3 |
| element halves the gas side in the denominator | RED x2 |
| element takes the jet ratios against the COOLANT density | RED |
| Fig. 10's A/B and A/C curves swapped | RED (Table 3) |
| the skew correction IS applied to the committed data | RED (Table 3) |

**Two results worth recording, because the first read was wrong both
times.**

The heat-up perturbation first looked GREEN, because it was run only
against the acceptance test, whose bands are deliberately wide -- plate C
moves to +0.1% and plate B to +19.5%, both still inside. Run against the
whole file it is red three ways, including the element/runner agreement
test. The lesson is about the falsification, not the code: a wide
acceptance band is not the binding test, and checking only it understates
coverage.

The skew perturbation IS caught, but only INCIDENTALLY. Correcting the
two scored curves breaks Table 3 because the other six stay raw. A
UNIFORM correction of all eight would slip through: it moves eta by at
most 0.005, and both the acceptance bands and Table 3's ratios are wider
than that. So two exact coordinates per scored curve are now pinned in
`test_fig10_calibration_and_why_the_skew_is_not_corrected` as the literal
record of what was digitised.
