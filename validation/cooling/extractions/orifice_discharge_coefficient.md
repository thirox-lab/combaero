# Extraction: orifice discharge coefficient with corner radius, length and crossflow

**Status: CONFIRMED (2026-09-23). I1 decided in favour of the full chain,
including the crossflow term; I2 and I3 follow from it. Cleared for
implementation, to be validated against the source's own figures. I5 (Wu)
deferred and out of scope.**

Every equation below was checked against the printed page image at 400 dpi, not
the OCR text layer, which garbles all three sources (McGreehan and Schotsch's
opens with "77je discharge coefficient"; Rohde's runs two-column nomenclature
into nonsense). Rendered with `pdftoppm -r 400`.

Tracked by #375. Surfaced from #337 (item 24b), under #339.

## Sources, pinned

> McGreehan, W.F. and Schotsch, M.J. (1988). Flow Characteristics of Long
> Orifices With Rotation and Corner Radiusing. *ASME J. Turbomachinery*,
> **110**(2), 213-217. Paper No. 87-GT-162. DOI 10.1115/1.3262183.
> `sources/orifices/mcgreehan_schotsch_1988.pdf` (gitignored, copyrighted).
>
> **Note the title.** #375 and the earlier session record both call this
> "...With Rotation and Crossflow". The printed title is "**...With Rotation
> and Corner Radiusing**". Crossflow is covered, but under the paper's own name
> for it, *relative tangential velocity*.

> Rohde, J.E., Richards, H.T. and Metzger, G.W. (1969). Discharge Coefficients
> for Thick Plate Orifices With Approach Flow Perpendicular and Inclined to the
> Orifice Axis. NASA TN D-5467, Lewis Research Center, October 1969.
> `sources/orifices/Rohde_NASA_TN_D-5467.pdf` (gitignored).

> Wu, D., Burton, R. and Schoenau, G. (2002). An Empirical Discharge
> Coefficient Model for Orifice Flow. *International Journal of Fluid Power*,
> **3**(3), 13-18. `sources/orifices/wu_burton_schoenau_2002_orifces.pdf`
> (gitignored).

**What each one actually is**, which differs from what #375 assumed:

| source | provides | implementable as-is |
|---|---|---|
| McGreehan and Schotsch | a **closed-form composite `C_d` chain** over Reynolds number, inlet corner radius `r/d`, orifice length `L/d`, and relative tangential (crossflow) velocity `U_1/V_i` | **yes** -- this is the model |
| Rohde | **charts only**: `C_d` against a velocity head ratio, one curve per main-duct Mach number. An explicit chart-lookup design procedure, no correlation | **no** -- data and validity statements only |
| Wu et al. | a closed-form `C_d(Re)` valid to `Re -> 0`, with coefficients fitted **per orifice type** to hydraulic valve data | form yes, constants are for valves |

Later addition (2026-11-08): Hay, Lampard and Benmansour (1983), crossflow on
**both** sides of the hole at 30 and 90 deg -- data only, unscored, the one
measured set on a hole-angle term and on crossflow direction. See
[hay_1983_crossflow_cd.md](hay_1983_crossflow_cd.md).
Gritsch, Schulz and Wittig (1998) add cylindrical, fan-shaped and laidback
fan-shaped 30 deg holes with coolant and mainstream crossflow (to Ma_m 1.2),
Cd on the inlet (metering) area -- data only, unscored, the one measured set
on a shaped-hole Cd. See
[gritsch_1998_shaped_hole_cd.md](gritsch_1998_shaped_hole_cd.md).

---

## For the reviewer

Five items need judgement. Everything else in this document is a transcription
fact, checked and reproduced below.

| # | what to decide | what turns on it |
|---|---|---|
| **I1** | **Adopt McGreehan and Schotsch as combaero's jet-plate `C_D`, replacing Florschuetz's fixed `0.79`?** It is a composite of *other people's* fits assembled by GE for network flow analysis -- the Re baseline is Miller and Kneisel's, the `r/d`, `L/d` and `U_1/V_i` terms are fitted to Benedict, Cusick, Lichtarowitz, Rohde, Grimm, and Meyfarth and Shine. Accepting it means accepting a chain of secondary fits. | Large. A predicted `C_D` runs `0.66-0.84` across `t/d = 0.5-3`, against one fixed `0.79`. `C_D` enters Florschuetz's `beta` (item 21 of `han_impingement.md`), so it moves every `Gc/Gj`. **Evidence for: check H below** -- the chain independently lands on `0.79` to 0.6% at `t/d = 1`, and stays inside Florschuetz's measured `0.73-0.85` band across `t/d = 1-3` and `Re_j = 5e3-7e4`. That is a different author, a different rig and a different decade agreeing. |
| **I2** | **Does the `U_1/V_i` term apply to a jet plate's crossflow, or only to rotation?** The paper's own answer is yes: p.215 names "an orifice in the side of a duct" as a source of relative tangential velocity, and Fig. 4's data includes Rohde's static side-of-duct orifice, not only rotating rigs. But the Summary hedges for *angled* orifices -- "the present method of solution is to assume that the `C_d` is influenced the same way ... Additional work is required to check this assumption against available data." | This is the whole point of #375. If `U_1/V_i` does not transfer, we get a geometry-only `C_D` and no crossflow dependence, and Florschuetz's `Gc/Gj` keeps a constant `C_D`. |
| **I3** | **Keep the non-monotonic bump?** Eq. (17) puts `C_d` **above** its no-crossflow value at small crossflow -- peak `+5.7%` at `U_1/V_i = 0.085` for a `0.60` baseline. I assumed a transcription error and was wrong: it is drawn in Fig. 4 and the data show it (check F). It is physical -- Rohde's result 2 is that slanting the orifice into the flow *increases* `C_d`. | A local maximum in `C_d` against crossflow means a non-monotonic solver residual in exactly the coupling `Gc/Gj` already makes implicit. Worth knowing before it presents as a convergence oddity. My recommendation: keep it, pin it with a test, and note it where the solver can see it. |
| **I4** | **Accept that Rohde contributes no correlation?** It is charts plus a "read the curve" design procedure. Its value here is (a) the validity statements below, several of which directly license what #375 wants to do, and (b) digitisable validation data. | If the reviewer wants Rohde as a *model*, Figs. 5, 7 and 8 need digitising and a fit -- a much larger scope, and one McGreehan and Schotsch already did. |
| **I5** | **Is Wu's low-Re extension in scope now, or later?** McGreehan and Schotsch Eq. (8) is floored at `Re >= 10 000`. Wu's form runs to `Re -> 0`, but its published constants are fitted to Merritt's sharp-edged curve, a spool valve and a needle valve -- hydraulic oil, valve geometry. Using them for a gas-turbine jet plate extrapolates across application, not just across Reynolds number. | My recommendation: **defer**, and implement Wu separately and labelled when something needs it, never blended silently into the McGreehan and Schotsch chain. A jet plate below `Re = 10 000` is a real case, but it deserves its own evidence rather than a borrowed valve fit. |

### What I still need sent

Nothing to proceed on I1-I4. For I5, a jet-plate `C_D` measurement below
`Re = 10 000` would settle whether Wu's constants transfer; absent one, deferring
is the honest call.

---

## Extracted items -- McGreehan and Schotsch (1988)

Equation numbers are the paper's own.

| # | item | content | status |
|---|---|---|---|
| 1 | Eq. (2), basic orifice equation | `W = C_d Y A_a sqrt(2 g_c P_t1/(R T_t1) (P_t1 - P_s2))` | confirmed |
| 2 | Eq. (3), approach-velocity factor | `k = 1/sqrt(1 - (d/D)^4)`; "for flow in a pipe the normal velocity of approach correction". Not used for plenum-to-plenum (`beta = 0`) | confirmed |
| 3 | Eq. (4), orifice expansion factor | `Y_o = 1 - 0.41 (P_t1 - P_s2)/(gamma P_t1)` | confirmed |
| 4 | Eq. (5), nozzle expansion factor | `Y_n = [S^(2/gamma) (gamma/(gamma-1)) ((1 - S^((gamma-1)/gamma))/(1 - S))]^(1/2)`, `S = P_s2/P_t1` | confirmed |
| 5 | Eq. (6)/(7), expansion-factor blend | `Y = (1-X) Y_o + X Y_n` for `C_d > 0.82`, with `X = 8.333 (C_d - 0.82)`. `X` reaches 1 at `C_d = 0.94` | confirmed |
| 6 | Eq. (8), Reynolds baseline, orifice | `C_d:Re = 0.5885 + 372/Re`, plenum-to-plenum (`beta = 0`), **valid `Re >= 10 000`**. Attributed to Miller and Kneisel [8] | confirmed |
| 7 | Eq. (9), Reynolds baseline, nozzle | `C_d:Re = 0.9981 - 4.73/sqrt(Re)`, **valid `Re >= 10 000`** | confirmed |
| 8 | high-Re asymptotes | 0.5885 (orifice), 0.9961 (nozzle) as printed. **0.9961 disagrees with Eq. (9)'s own limit of 0.9981** -- see item 8a | confirmed as printed |
| 8a | the 0.9961/0.9981 discrepancy | Eq. (9) tends to 0.9981, not 0.9961, as `Re -> inf`. A 0.2% printing slip in the running text. Nothing depends on it: the chain uses Eq. (9), not the quoted asymptote | suspect (text), harmless |
| 9 | reference point | "A baseline `C_d` of 0.60 for a sharp-edged orifice is used as a reference point at `Re = 3.2 x 10^4`". The `0.6` divisors throughout Eqs. (1) and (17) are this number | confirmed, and reproduced -- check A |
| 10 | Eq. (10), Schoder and Dawson corner radius | `dC_d/C_d = 3.10 (r/d)`, "a good approximation for `0 < r/d < 0.1`". Cited [13], **not** the method used | confirmed |
| 11 | Eq. (11), corner-radius correction | `C_d:r = 1 - f (1 - C_d:Re)` | confirmed |
| 12 | Eq. (12), corner-radius function | `f = 0.008 + 0.992 e^(-5.5 (r/d) - 3.5 (r/d)^2)` | confirmed |
| 13 | `r/d` data range and limit | Cusick [3] covers `0.1 < r/d < 0.36` for quadrant-edge orifices; Benedict et al. [14] `r/d < 0.0035`. An ASME nozzle is the limiting `C_d`, reached at `r/d = 0.82`; beyond that "no further benefit for inlet radiusing is achieved" | confirmed |
| 14 | Eq. (13), length correction | `C_d:r,L = 1 - g (1 - C_d:r)` | confirmed |
| 15 | Eq. (14), length function | `g = [1 + 1.3 e^(-1.606 (L/d)^2)](0.435 + 0.021 L/d)` | confirmed |
| 16 | Eq. (15), combined `r/d` and `L/d` | `(C_d:r)' = 1 - g (1 - C_d:r)`, with **`g` evaluated at `r/d`** in place of `L/d`. Rationale given: the corner radius reduces inlet separation and so reduces the long orifice's dynamic-pressure-recovery benefit | confirmed |
| 16a | **the prose and the formula of Eq. (15) point in opposite directions** | The rationale given is that the corner radius "decreases the dynamic pressure recovery benefit of the long orifice" -- i.e. less length benefit. The printed formula does the reverse. Eqs. (13)+(15)+(16) collapse to `Cd = 1 - g(L/d - r/d) g(r/d) (1 - Cd:r)`, and `g(L/d - r/d) g(r/d) < g(L/d)` for every `r/d > 0`, so the combined correction yields *more* length benefit, not less (at `L/d = 2`: `0.478` falls to `0.447` at `r/d = 0.3` and `0.402` at `r/d = 0.5`). Net `Cd` is still monotone in `r/d` and still approaches the nozzle limit, so nothing downstream misbehaves | confirmed as printed; the **formula is implemented**, the prose is not. Flagged so a future reader does not "fix" the code to match the rationale |
| 17 | Eq. (16), effective length | **Printed `(L/d)' = L/D - r/d`.** Read as `L/d - r/d`; the capital `D` is a typo -- see item 17a | suspect as printed, resolved by derivation |
| 17a | why `L/D` cannot be meant | Three independent reasons. (i) The prose says "the inlet radius should be subtracted from the total length", i.e. both terms in units of `d`. (ii) Nomenclature defines `D` as *pipe diameter*; this chain is plenum-to-plenum, where `beta = 0` and `D` is undefined, giving `(L/d)' = -r/d < 0`. (iii) `(L/d)'` must reduce to `L/d` at `r/d = 0`, which only `L/d - r/d` does. Glyph checked at 400 dpi: the printed character is a capital `D`, distinct from the `d` in `r/d` on the same line, so this is the paper's error and not a misread | resolved |
| 18 | Eq. (17), relative tangential velocity | `C_d:r,L,U = C_d:r,L (C_1 + C_2 C_3)` | confirmed |
| 19 | Eq. (17) terms | `C_1 = e^(-R_v^1.2)`; `C_2 = 0.5 R_v^0.6 (C_d:r,L/0.6)^(-0.5)`; `C_3 = e^(-0.5 R_v^0.9)` | confirmed, and the `-0.5` sign falsified -- check G |
| 20 | velocity ratio parameter | `R_v = (U_1/V_i)(C_d:r,L/0.6)^(-3)` | confirmed |
| 21 | `U_1/V_i` sources | Rotating orifice with stationary inlet air, **or an orifice in the side of a duct** (crossflow). Data compiled from Rohde [10], Grimm [11], Meyfarth and Shine [12] | confirmed |
| 22 | the crossflow caution | "Caution must be taken not to use the total pressure based on the tangential velocity when calculating ideal orifice velocity with inlet crossflow (orifice in the side of a duct). This will significantly alter the resulting `C_d`." `V_i` is built from **static** pressure and temperature consistent with flow through the orifice | confirmed -- a correctness trap for the implementation |
| 23 | Rohde adjusted | "Rohde's results have been adjusted using static parameters", Rohde having calculated orifice velocity from duct *total* conditions | confirmed |
| 24 | Eq. (1), the superseded method | `C_d = (C_d:r/0.60)(C_d:r,L/0.60) x 0.60`. Explicitly the **original** method, which "can result in discharge coefficients greater than 1.0, which are illogical". Superseded by Eqs. (8)-(17) | confirmed -- **must not be implemented** |
| 25 | accuracy claim vs Rohde | "The results of Rohde generally show trends that are consistent ... although his basic values are lower for `r/d` and `L/d > 0`." No error statistic is given anywhere in the paper | confirmed -- no stated tolerance exists |
| 26 | Eq. (18), flow function | `W sqrt(T_t1)/P_t1 = C_d Y A_a sqrt(2 g_c/R (1 - P_s2/P_t1))` | confirmed |

## Extracted items -- Rohde et al. (1969)

| # | item | content | status |
|---|---|---|---|
| 30 | correlating parameter | "velocity head ratio" `(P_T - p_j)/(P_T - p_d)`: velocity head of the orifice jet over velocity head in the main duct. `p_j` = jet exit static, `p_d` = main duct static, `P_T` = main duct total | confirmed |
| 31 | form of the result | Charts of `C_d` against velocity head ratio, one curve per main-duct Mach number. **No correlation equation anywhere in the report** | confirmed |
| 32 | tested ranges | `t/d = 0.51` to `4.00`; `d = 0.059-0.128 in` (0.150-0.325 cm); `t = 0.06-0.25 in`; inlet edge sharp to `r = 0.030 in`; main duct `M = 0` to `0.65`; `p = 20-80 psia`; `T` ambient to ~1000 F (811 K); orifice axes at 45 deg and 90 deg to the duct | confirmed |
| 33 | **multiple-orifice interference** | Summary of Results 4: "The longitudinal interference of multiple orifices is **not significant above an orifice center distance to orifice diameter ratio of 1.5** for an orifice with a wall thickness to orifice diameter ratio of 1." Abstract and Summary both list multiple-orifice interference among the effects found **negligible** | confirmed -- **this licenses using a single-hole `C_d` for a jet array** |
| 34 | predominant factors | Summary of Results 1: approach Mach number, static pressure differential, `t/d`, and inlet edge radius | confirmed |
| 35 | inclined orifices | Summary of Results 2: "Slanting the orifice axis in the direction of flow **significantly increases** the discharge coefficient" | confirmed -- independent support for item I3's bump |
| 36 | negligible effects | Summary of Results 5: air temperature, pressure level, entrance length, orifice surface finish. Result 3: duct wall curvature is a small effect | confirmed |
| 37 | extrapolation caveat | Appendix: "The variation of placing the orifice on an angle or the inlet edge condition have such a marked effect on discharge coefficient, that its extrapolation to `t/d` ratios other than 1.0 would be questionable" | confirmed -- bounds how far item 35 can be pushed |
| 38 | why the report exists | Introduction: predicting flow in internally cooled turbine blades and vanes; thick plate orifices "employed in the impingement cooling concept". Same application as #337 | confirmed |

## Extracted items -- Wu, Burton and Schoenau (2002)

| # | item | content | status |
|---|---|---|---|
| 40 | Eq. (3), two-parameter form | `C_d = C_dinf (1 - e^(-(delta/C_dinf) sqrt(Re)))` | confirmed |
| 41 | Eq. (4), general form | `C_d = C_dinf (1 + a e^(-(delta_1/C_dinf) sqrt(Re)) + b e^(-(delta_2/C_dinf) sqrt(Re)))` | confirmed |
| 42 | sharp-edged orifice, Merritt's curve | Table 1 inputs `C_dinf = 0.61`, `delta = 0.23`, `C_dm = 0.69`, `sqrt(Re_m) = 11`; outputs `a = 1.07`, `b = -2.07`, `delta_1 = 0.077`, `delta_2 = 0.15` | confirmed, and internally cross-checked -- check J |
| 43 | Eq. (15), the fitted result | `C_d = 0.61 (1 + 1.07 e^(-0.126 sqrt(Re)) - 2.07 e^(-0.246 sqrt(Re)))` | confirmed |
| 44 | spool orifice fit (Fig. 4) | `C_d = 0.63 - 0.625 e^(-0.21 sqrt(Re)) + 0.005 e^(-5.0 sqrt(Re))` | confirmed; the `-5.0` exponent is the least legible digit on the page |
| 45 | needle valve fit (Fig. 5) | `C_d = 0.75 - 2.47 e^(-0.22 sqrt(Re)) + 1.72 e^(-0.28 sqrt(Re))` | confirmed |
| 46 | `delta` for a sharp edge | "For a sharp edged orifice, `delta` is approximately 0.2" | confirmed |
| 47 | what the data are | Hydraulic fluid power: Merritt's textbook curve, a Brand Hydraulics EFC12-10-12 fixed orifice, a Deltrol EN-35 needle valve. **No gas-turbine or jet-plate data** | confirmed -- the basis for I5 |

---

## Checks that could have failed

Run by `tmp/ms_checks.py` (transcribe-then-verify; the script is scratch, and
the checks it encodes move into the test suite on implementation). Each has an
anchor outside the equation being tested.

| check | anchor | result |
|---|---|---|
| **A** | Eq. (8) evaluated at `Re = 3.2e4` against the "baseline 0.60" stated in *prose on a different page* | `0.60013` vs `0.60` -- **0.02%** |
| **B** | Eq. (12) against Schoder and Dawson's Eq. (10), *a different author's correlation* | 0.73% at `r/d = 0.02` and `0.05`, 1.31% at `r/d = 0.10` (Eq. 10's stated upper limit) |
| **C** | Eq. (11)+(12) extrapolated to `r/d = 0.82` against Eq. (9)'s nozzle asymptote -- *a different equation, fitted to different data* | `0.99628` vs `0.9981` -- **0.18%**. The corner-radius fit, pushed to the radius at which the paper says a nozzle is reached, lands on the paper's own nozzle correlation |
| **D** | identities the fitted constants must satisfy: `g(0) = 1`, `f(0) = 1`, Eq. (17) `-> C_d:r,L` at `U_1/V_i = 0` | 0.05%, exact, exact |
| **E** | Eq. (13)/(14) against the drawn curve of Fig. 3 | 0.0-1.1% for `L/d >= 1`; ~3% at `L/d = 0.5`, on the steepest part of the curve where a 0.1 error in reading `L/d` moves `C_d` by 0.017 |
| **F** | Eq. (17) against the drawn curves of Fig. 7, and the bump against Fig. 4 at 400 dpi | 0.1-0.7% at `U_1/V_i <= 2`, ~3% at `U_1/V_i = 4`. **Bump: computed peak `0.6344` at `U_1/V_i = 0.085`; drawn peak ~`0.635` at ~`0.1`** |
| **G** | falsification of the `C_2` exponent sign: substituting `(C_d:r,L/0.6)^(+0.5)` for `^(-0.5)` | gives `C_d = 1.0025` against a `C_d:r,L = 1.0` baseline -- a correction factor exceeding 1, which the paper explicitly calls out as the defect of the superseded Eq. (1). The sign is confirmed by the failure of its alternative |
| **H** | **the #375 anchor.** The chain at `r/d = 0`, `t/d = 1`, `Re_j = 1e4` against Florschuetz et al. (1981)'s recommended jet-plate `C_D = 0.79` -- *a different paper, rig and decade, and the number #337 currently hard-codes* | `0.7848` vs `0.79` -- **0.6%**. Across `t/d = 1-3` and `Re_j = 5e3-7e4` the chain gives `0.766-0.839`, inside Florschuetz's measured `0.73-0.85` band (Table 1). At `t/d = 0.5` it gives `0.66-0.72`, below the band |
| **J** | Wu Table 1's `delta_1`, `delta_2`, `C_dinf` against the exponents printed in Eq. (15), transcribed separately | `0.077/0.61 = 0.1262` vs `0.126`; `0.15/0.61 = 0.2459` vs `0.246` |

**What check H does not show.** It is an agreement at one geometry, not a
validation. `t/d = 0.5` falls outside Florschuetz's band, and Florschuetz's own
`0.73-0.85` spread is wider than the difference between a predicted `C_D` and a
fixed `0.79` over much of the range. H is strong evidence that the chain is
*calibrated for this configuration*; it is not evidence that a varying `C_D`
predicts any measured jet-array flow split better than `0.79` does. Only
scoring `Gc/Gj` against data can show that, and no such data is in hand.

---

## Modelling decisions

**D1: implement the chain as separate stages, not one composite function.**
`orifice.h` already exposes per-correlation free functions over
`OrificeGeometry`/`OrificeState`. Each of Eqs. (8), (11)-(12), (13)-(16),
(17) is independently citable and independently testable, and checks B, C, E
and F each bind to one stage. A single opaque `Cd_mcgreehan_schotsch()` would
make every one of those a test of the whole chain.

**D2: do not implement Eq. (1).** The authors supersede it in the same paper
and name its failure mode (`C_d > 1`). Recorded only so a later reader does not
find it in the paper and think it was missed.

**D3: `(L/d)' = L/d - r/d`, per item 17a**, with the paper's printed `L/D`
recorded in a comment as a known erratum. A reader who checks the code against
the paper will otherwise find a discrepancy and "fix" it.

**D4: keep Eq. (17) non-monotonic** (item I3), subject to the reviewer. It is
in the source, in the source's data, and independently supported by Rohde's
result 2. Suppressing it would be tuning the model to suit the solver.

**D5: `V_i` must be built from static conditions** (item 22). This is the one
place in the chain where a plausible implementation is silently wrong: using
duct total pressure inflates `V_i`, deflates `U_1/V_i`, and biases `C_d` high
by exactly the amount the paper warns about. Worth an explicit named argument
rather than a comment.

**D6: solver-facing, so C++ with an analytic `(f, J)`**, per the project's
`(f, J)` rule -- `C_D` feeds Florschuetz's `beta` and hence `Gc/Gj`, which is
inside the residual. The chain is a composition of elementary functions, so
forward-mode `DualN` applies directly with no new derivation; the established
recipe is in the C++ Jacobian port pattern. Note the derivative of Eq. (17) is
where `R_v^0.6` and `R_v^0.9` produce an infinite slope at `R_v = 0`, so the
`U_1/V_i -> 0` limit needs a guard that a finite-difference check will not
reveal on its own.

**D8: `U_1/V_i` is a SUPPLY-side ratio and must never be fed Florschuetz's
`Gc/Gj`.** Item 21 and item 22 both define `U_1` as the *inlet* relative
tangential velocity -- the approach flow upstream of the orifice, on the plenum
side of a jet plate. Florschuetz's `Gc/Gj` is spent-air crossflow accumulating
in the impingement channel, *downstream* of the holes. They are opposite faces
of the same plate and are not interchangeable:

- **Florschuetz's plate is plenum-fed**, so its approach velocity is
  negligible and the correct evaluation is `U_1/V_i = 0`. That is precisely
  the configuration check H anchors at 0.6%, so the anchor is valid.
- **A duct-fed orifice** -- Rohde's configuration, an orifice in the side of a
  passage, and a real combaero case -- is where `U_1/V_i > 0` applies.

Both belong in the implementation; the guard is that the crossflow argument is
named for what it is (approach/supply-side) and documented against this trap,
rather than being a bare `crossflow_ratio` that a caller will eventually feed
`Gc/Gj`. Feeding it `Gc/Gj` would be wrong in a way nothing would catch: it
produces a plausible number, biased low, from the wrong face of the plate.

**D9: keep an explicit constant-`C_D` path.** Two exist and both need a small
fix rather than a new mechanism:

- Python, jet-plate side: `ImpingementModel.C_D` is already a plain float
  defaulting to `FLORSCHUETZ_1981_DEFAULT_CD`. It stays a float, and the
  correlation is offered as a function a caller evaluates and passes in. No
  union type and no auto-selection, so today's behaviour and default are
  unchanged and the switch is visible at the call site.
- C++, orifice side: `ConstantCdCorrelation` exists and takes a value, but
  `make_correlation(CdCorrelation::Constant)` hardcodes its `0.61` default with
  no way to pass one -- the enum can be selected but the constant cannot be
  set. That is a real gap in the switch, not a style point.

**D10: `Y` belongs to the INCOMPRESSIBLE path only.** Settled by reading the
code, not by preference. combaero's `regime='compressible'` calls
`combaero::nozzle_flow` with `A_eff = Cd E A` (`src/solver_interface.cpp`), a
full isentropic solve with a choked branch and real-gas properties. Eq. (5)'s
`Y_n` **reproduces that solve to 0.008%** over `0.6 <= S <= 0.99` -- measured,
not assumed. So `nozzle_flow` already *is* the expansion factor, and applying
`Y` on top of it corrects for compressibility twice (worth >8% at `S = 0.7`,
pinned by `test_applying_y_on_top_of_nozzle_flow_would_double_count`).

`Y` therefore does exactly one job: it extends the incompressible form
`m_dot = Cd A sqrt(2 rho dP)` to finite pressure ratio. `orifice.h`'s module
comment, which claimed the incompressible correlations combine with
`compressible.h` without saying who owns the expansion term, is updated.

That 0.008% agreement is also the strongest single check on the Eq. (5)
transcription: a closed form from a 1988 paper against combaero's own
independent numerical solve.

**D11: saturate Eq. (7)'s blend weight smoothly, not with a hard clamp.**
The paper prints `X = 8.333 (Cd - 0.82)` as a bare ramp: it says the blend
applies "for `C_d > 0.82`", `X` reaches 1 at `C_d = 0.94`, and no upper bound
is given. Unclamped, our own chain reaches `X = 1.45` at `Cd = 0.994`, giving
`Y = -0.45 Y_o + 1.45 Y_n` -- an extrapolation past a nozzle with no backing.
So a clamp is needed, and **the clamp is a decision the source did not make**.

A HARD clamp was rejected on measurement. It puts `dY/dCd` at *exactly zero*
outside `[0.82, 0.94]` -- and `Cd > 0.94` is a real design point (`r/d >= 0.2`
at `t/d = 2`) -- plus a discontinuous jump of **0.734** in `dY/dCd` at both
knees. The discontinuity is the disqualifying part: it breaks the
finite-difference agreement any `(f, J)` pair has to satisfy, which is
squarely in the way of #383.

The chosen form is a smooth saturation whose `eps -> 0` limit is *exactly* the
paper's hard clamp, so the fidelity case is available rather than lost:

| | hard clamp | smoothed, `eps = 0.05` |
|---|---|---|
| `dY/dCd` outside the blend | exactly 0 | non-zero everywhere reachable |
| change in `dY/dCd` across a knee, sampled at `h` | **0.73353 for every `h`** | halves each time `h` halves |
| worst `|Y_soft - Y_hard|` | -- | **0.0042, i.e. 0.47% of a typical `Y`** |

Continuity is tested by that scaling rather than against a threshold: a
genuine discontinuity holds the same jump however closely it is sampled, a
continuous derivative's shrinks with the interval.

**`eps = 0.05` is chosen by measurement, not judgement.** Two independent
properties bound it from opposite directions, and the window between them is
narrow:

| `eps` | width in `Cd` | cost in `Y` | continuity ratio | min `dY/dCd` | verdict |
|---|---|---|---|---|---|
| 0 (hard) | 0 | 0 | 1.00 | **0** | floor + kink |
| 0.005 | 0.0006 | 0.03% | 0.90 | 1.4e-6 | effectively floored |
| 0.02 | 0.0024 | 0.14% | 0.60 | 2.2e-5 | continuity fails |
| **0.05** | **0.0060** | **0.35%** | **0.52** | **1.4e-4** | **ok** |
| 0.10 | 0.0120 | 0.68% | 0.51 | 5.4e-4 | fidelity fails |
| 0.30 | 0.0360 | 1.95% | 0.50 | 4.8e-3 | fidelity fails |

Bisected walls: continuity holds for `eps >= 0.0304`, fidelity for
`eps <= 0.0976`. The admissible window is **[0.0304, 0.0976]** and the default
sits 29% across it. Falsified in both directions -- setting the constant to
0.005 or 0.02 turns the dead-zone and window tests red, 0.15 turns the
fidelity and window tests red, and the failing test names which wall was hit.

A note on the lower wall, because it is easy to get wrong: "non-zero
derivative" is not the same as "usable derivative". The first version of the
dead-zone assertion tested `!= 0`, which passes at `eps = 0.005` where the
smallest slope is `1.4e-6`. It now requires room (`> 5e-5`), which makes it a
second, independent lower wall agreeing with the continuity one.

**D14: Cd's own solver treatment, and the two hazards it had to clear.**
Analytic `(f, J)` for the chain landed with `cd_and_derivatives`, returning
`(Cd, dCd/dRe, dCd/d(U1/Vi))`. `r/d` and `L/d` are geometry and carry no
partials. Two hazards, deliberately treated DIFFERENTLY because measurement
said they are not the same kind of problem.

**Hazard 1, the crossflow singularity -- regularised in the MODEL.**
Eq. (17)'s `Rv^0.6` and `Rv^0.9` have unbounded slope at `Rv = 0`, and
`U1/Vi = 0` is both the default and the physically right value for a
plenum-fed jet plate, so the singularity sits where the solver lives. The
input is regularised, `U1/Vi -> sqrt((U1/Vi)^2 + eps^2)`.

`rv_smooth_eps = 3e-5` is bounded by two walls, both sourced:

| wall | criterion | value |
|---|---|---|
| lower | worst `|dCd/du|` within 10x the physical scale (0.663, measured on the exact chain over `u >= 0.01`, the range Figs. 4-6 carry data for) | `eps >= 1.86e-5` |
| upper | departure from Eq. (17) within 10% of the scatter the correlation sits in (`+/-0.02` in `Cd`; **the paper states no error statistic at all** -- item 25) | `eps <= 4.15e-4` |

The default sits **15% across**, below where `x_smooth_eps` sits in its own
window (29%), deliberately. Reviewer's operational experience is that
OVER-smoothing is the worse failure on this solver -- a residual that no
longer matches the physics stalls the solve, and that bit harder than a stiff
Jacobian ever did. 3e-5 keeps 5x headroom against the fidelity wall while
holding margin on the stiffness one (8.2x of 10x); 2e-5 smooths marginally
less but sits at 9.7x with no room for a re-measurement to move it.

**Two sampling errors were made reaching that number, and both flattered the
answer.** The first linear sweep stepped `u` by 0.002 and never sampled below
it -- the whole singularity lives under 1e-3 -- and reported a "plateau at
1.451" that does not exist; the true relation is a clean `eps^-0.4` power law,
exactly as `Rv^-0.4` predicts, so every `eps` trades derivative against
fidelity and there is no free choice. The second measured the fidelity cost
over `u >= 0.01` only, excluding `u = 0` -- the default, and where the
regularisation bites hardest -- reporting 3.3e-7 where the truth is 8.6e-4,
three orders out, which moved the upper wall by a factor of 21. **A log grid
including the endpoint is not optional here.** Both were caught by a test, not
by re-reading.

**Hazard 2, the Re dead zone -- treated in the JACOBIAN ONLY.**
Below `re_min` the value is held at the floor, which is right: Eq. (8) is
outside its validity and divergent there. That leaves `dCd/dRe` exactly zero
over a third of a typical scan range, plus a kink at the boundary.

Softening the VALUE was measured and rejected: a soft-max floor recovers only
16% of the live slope at `eps = 1e4` in Re units while moving `Cd` by 7.1e-3,
and it would distort a region where the correlation does not apply. Instead
the derivative alone is continued from `re_min` -- implemented by seeding
`max(Re, re_min)` as the dual unknown, so the value stays floored while the
partial comes out as the live slope at the floor. No reported `Cd` changes,
and `dCd/dRe` becomes continuous across `re_min`.

The asymmetry is the finding: hazard 1 needed the model regularised, hazard 2
needed only the Jacobian. Treating them alike would either distort valid
physics or leave the dead zone in place.

FD cross-checked to 3.7e-7 relative on `dCd/dRe` and ~1e-9 on `dCd/du`.
`u = 0` is excluded from the FD comparison and separately reasoned: the
`max(U1_over_Vi, 0)` clamp makes a central difference one-sided there, so FD
measures a forward slope that cannot be compared; the analytic zero is correct
because `sqrt(u^2 + eps^2)` has zero slope at the origin.

**The walls hold across the geometry range, and that was checked rather than
assumed.** Both were bisected at one operating point, and there was a
specific reason to doubt them elsewhere: `Rv` carries `(cd_base/0.6)^-3`, so
a rounded long hole at `Cd ~ 0.95` sees the same `U1/Vi` as a 4.5x smaller
`Rv`, and the regularisation is applied to `U1/Vi` rather than `Rv`. Swept
over `Re` 1e4-1e6, `r/d` 0-0.2, `L/d` 0.5-5 (`cd_base` 0.76-0.95):

| | across the sweep | wall |
|---|---|---|
| ratio to physical scale | 7.8 - 8.3x | 10x |
| cost | 3.2e-4 - 4.2e-4 | 2.0e-3 |

The ratio is nearly geometry-INVARIANT, and the reason matters: the `Rv`
stretching moves the regularised peak and the physical maximum together, so
they cancel. That is what makes a single-point bisection a safe way to set
`eps` here -- and the prediction that high `cd_base` would be stiffer was
wrong in direction (it is 7.8x there against 8.2x at `cd_base` 0.77). Pinned
by `CrossflowSmoothingWallsHoldAcrossTheGeometryRange`, which also asserts
the invariance itself, since if the ratio stops being flat the single-point
method stops being valid.

Falsified four ways: `eps` below the lower wall, `eps` above the upper wall,
and each derivative partial dropped in turn.

**D13: every bound in the expansion factor is soft, and that is scanned for,
not assumed.** A sweep for zero-derivative regions and derivative
discontinuities across `Cd`, `S` and `Re` caught a bound this document had
already claimed was handled: an earlier revision smoothed the choke point but
left a hard `min(S, 1)`, putting a kink of `~0.455` in `dY/dS` at `S = 1` and
a floor above it. `S > 1` is reverse flow and a solver iterate reaches it. The
scan found it; reading the code had not. Pinned by
`PressureRatioIsSmoothAtBothBounds`, falsified against the original bug.

That scan is now `validation/solver_smoothness.py`, reusable and tested
against hazards with known analytic answers; `test_orifice_expansion.py` runs
it over the whole `Y` surface as a standing guard, and over the `Cd` chain
with #383's two known hazards named in `allow=` so the rest stays guarded and
a NEW hazard still fails. Three defects in the detector itself were found and
fixed by those reference cases: it measured floor width linearly on log scans
(hiding a flat decade), it skipped domain-edge points (hiding the one hazard
that matters most, since `U1/Vi = 0` is both an edge and the default operating
point), and its first boundary test confused "derivative still changing" with
"derivative unbounded", flagging `log(x)`.

The same scan re-confirms the two hazards already tracked in #383 and finds no
others: `dCd/dRe` is exactly zero below `Re = 1e4`, and `dCd/d(U1/Vi)` at
`U1/Vi = 0` is not merely discontinuous but divergent -- its jump *grows* from
24.7 to 43.1 as the sampling interval is quartered, where a kink would hold
constant and a smooth point would fall by four.

**D12: saturate `S` at the critical pressure ratio, for the same reason.**
Below `S* = 0.5283` (for `gamma = 1.4`) the isentropic form turns over and
predicts *decreasing* flow -- the classic non-physical branch, and the proxy
mass flux `Y_n sqrt(1-S)` was confirmed to peak exactly at `S*`. Y saturates
there. The paper gives no choked branch (Eq. (5) is a subsonic isentropic
form), so this closure is ours; it is the textbook one, and it is smoothed
with the same mechanism so it does not reintroduce a dead `dY/dS`.

**D7: defer Wu** (item I5). Implement only when a case needs `Re < 10 000`, and
then as its own labelled correlation with its valve provenance stated, never
blended into this chain.

---

## The Rohde conversion

Rohde plots `Cd` against a velocity head ratio `VHR = (P_T - p_j)/(P_T - p_d)`
built from duct TOTAL conditions; Eq. (17) is written in `U1/Vi` with STATIC
parameters, and the paper states Rohde's results "have been adjusted using
static parameters" (item 23). Deriving that adjustment, with `p_d` the duct
static and `p_j` the jet exit static (both confirmed from the model drawing,
Fig. 4 of Rohde):

- `P_T - p_d = 0.5 rho U_1^2`, the approach velocity head
- `p_d - p_j = 0.5 rho V_i^2`, the static-referenced ideal jet head
- hence `U1/Vi = 1/sqrt(VHR - 1)`
- and, since `Cd = W/(rho A V_ideal)` with Rohde's `V_ideal` referenced to
  `P_T` instead of `p_d`, `Cd_static = Cd_Rohde * sqrt(VHR/(VHR - 1))`

**This derivation is recorded but deliberately NOT used to produce the
committed data.** It is divergent as `VHR -> 1` -- the `Cd` factor is 4.1x at
`U1/Vi = 4` -- so it is least trustworthy exactly where the crossflow term is
most interesting. Digitising the paper's Fig. 6 instead takes the conversion
from the authors who defined both coordinate systems. The derivation's role is
to make their replotting checkable, and to allow Rohde's other figures to be
brought in later if the reviewer wants them.

## What Figs. 5 and 6 actually test

Found while wiring the validation, and it changes what the comparison means:
**the curves in Figs. 5 and 6 are not chain predictions.** The paper anchors
them to "a set baseline point at `U1/Vi = 0`" (p.216) taken from Rohde's own
measurement. Those baselines are `0.64`, `0.73` and `0.88`; the chain predicts
`0.669`, `0.792` and `0.990` for the same geometries. The 0.73 label is not
even reachable by the chain -- it would need `Cd:Re = 0.48`, below Eq. (8)'s
0.5885 floor -- which is what proves the labels are measurements, not
predictions.

So Figs. 5 and 6 validate **Eq. (17) in isolation**, and the baseline chain is
validated separately by Figs. 2 and 3. Scoring the full chain against this
data measures the sum of two disagreements and blames both on crossflow. The
runner therefore has two modes, and `eq17` is the one that reproduces the
source's own validation.

## Measured agreement

Digitised: `validation/cooling/data/mcgreehan_schotsch1988/fig6_rohde_td0p51.csv`,
11 points, Fig. 6's open-square series (`(r/d, L/d) = (0, 0.51)`, baseline
0.64). Scored by `validation/cooling/orifice_runner.py`, pinned by
`python/tests/test_orifice_validation.py`.

| what | N | bias | RMS |
|---|---|---|---|
| **Eq. (17) at the figure's own 0.64 baseline** | 11 | **+5.95%** | **7.54%** |
| full chain, baseline not anchored | 11 | +13.32% | 15.45% |

The correlation reads **high** against Rohde, and increasingly so with
crossflow: `+3%` below `U1/Vi = 0.4`, `+16%` at `1.41`. That direction is the
paper's own statement that Rohde's "basic values are lower" (item 25), now
with numbers on it:

| `(r/d, L/d)` | Rohde's measured basic `Cd` | chain at `Re = 3.2e4` | gap |
|---|---|---|---|
| `(0, 0.51)` | 0.64 | 0.6692 | **+4.6%** |
| `(0, 4.0)` | 0.73 | 0.7925 | **+8.6%** |
| `(0.49, 1.06)` | 0.88 | 0.9902 | **+12.5%** |

The gap grows with `r/d` and `L/d`, exactly as item 25 says it should.

**The metric was falsified**, not just the change: perturbing `U1/Vi` moves RMS
to 26.4% at `x0.5` and 17.9% at `x2.0`, so the score is measuring the crossflow
term rather than an imposed quantity.

**A `x1.25` scaling of `U1/Vi` would cut RMS to 1.88%. It is not applied.**
Eq. (17) was fitted to Rohde, Grimm, and Meyfarth and Shine pooled; refitting
it to the single series it is being judged against is tuning a constant against
its own score. The `+5.95%` bias is the honest number and is what the test
band records.

## Harness gates the series had to pass

The dataset's own quality gates caught two defects in this digitisation that
the scoring did not, and both are worth recording because neither would have
shown up as a wrong number:

- **No verification card.** `verify.py` refuses to score an unchecked series.
  The card now records the axis ticks, the multiplier, and the point count
  read off the page *independently of where the detector put anything* -- 11
  squares counted on the render, 11 in the CSV.
- **The redrawn plot came out on a log ordinate.** Fig. 6's axes are linear.
  `test_axis_specs_are_declared_for_plotting` requires `x_axis_type` and
  `y_axis_type` to be declared rather than defaulted, for exactly this reason:
  a log redraw of a linear figure looks entirely plausible. Caught by the
  gate, not by inspection.

All 351 `verify.py` checks pass with the series included.

## Rohde Figure 10: the r/d term against data

Digitised by the reviewer (`validation/cooling/data/rohde1969/`): three inlet
edge conditions at `t/d = 1.06`, model 8, duct Mach 0.14. `r/d` from Table I,
which gives model 8 `d = 0.0615 in`, so the 0.012 in and 0.030 in radii are
`r/d = 0.195` and `0.488`. This is the only data in the project that tests
Eq. (12).

**Digitisation quality.** Counts were taken by two channels independently, the
automated count sealed to a file before the reviewer's data arrived:

| series | reviewer | sealed | |
|---|---|---|---|
| circles `r/d = 0.488` | 8 | 8 | agree |
| squares `r/d = 0.195` | 10 | 10 | agree |
| triangles `r/d = 0` | 12 | 12 | agree |

Agreement includes both series the sealed note had flagged as its
lowest-confidence calls. Calibration: x ticks within **0.59%**, y ticks within
**0.00185** in `Cd`. Ordering invariant (`circle > square > triangle` at every
velocity head ratio): **0 violations** in 9 checks.

**A y-skew, measured not guessed.** The corners file shows the canvas reading
high on its right-hand side -- mean `dy` `+0.00082` at the left edge,
`+0.02132` at the right, a skew of `+0.0205` in `Cd` across the log abscissa.
It is page rotation, not picking error: the y ticks, all picked at the left
edge, are good to `0.00185`. **Confirmed by a third channel, the page's own
words**: Rohde writes that the largest inlet radius reaches "as high as 0.94".
The circle series' raw maximum is `0.9569` (+1.8% against that text);
skew-removed it is `0.9362` (-0.4%). The raw picks are committed uncorrected
and the runner derives the correction from the committed corners file, so it
is reproducible rather than a typed constant.

### The finding: the r/d term does not capture its interaction with crossflow

The correlation reads **high** against Rohde, and the error grows sharply with
crossflow -- but not equally across the three series:

| series | VHR 1-10 | VHR 10-25 | VHR 25-60 |
|---|---|---|---|
| circles `r/d = 0.488` | +15.4% | +7.8% | +8.1% |
| squares `r/d = 0.195` | +27.2% | +6.4% | +4.0% |
| triangles `r/d = 0` | **+40.6%** | **+25.9%** | +8.9% |

**That ordering rules out the conversion as the cause.** The total-to-static
factor `sqrt(VHR/(VHR-1))` is identical for all three series at a given
velocity head ratio, so it cannot produce an error that is monotone in `r/d`.
The differential part of this disagreement is physics.

Equivalently, as a ratio test that is conversion-free by construction (at
fixed VHR the conversion maps all three to the same `U1/Vi` and scales all
three `Cd` identically, so the ratio between series cancels it):

| | measured `Cd(0.488)/Cd(0)` | predicted | |
|---|---|---|---|
| VHR 4 (`U1/Vi` 0.58) | 1.72 | 1.32 | **-23.5%** |
| VHR 45 (`U1/Vi` 0.15) | 1.22 | 1.26 | +3.3% |

Overall bias `-13.1%`, RMS `15.7%`. **Rounding the inlet protects an orifice
against crossflow degradation far more than Eq. (12) plus Eq. (17) predicts**,
and at low crossflow the ratio test agrees to `+3.3%`, so Eq. (12) itself is
sound -- it is the *interaction* that is missing. Physically reasonable:
rounding helps most when the flow has to turn hardest.

**Two independent confirmations that this is the model and not the data:**

1. At `VHR >= 25`, where the conversion factor is 1.013 and crossflow is weak,
   all three series converge on `+4%` to `+9%` -- matching the `+5.95%`
   measured from the wholly independent Fig. 6 digitisation, a different
   figure at a different `t/d` with no conversion used at all.
2. The disagreement is monotone in `r/d`, which no common factor can cause.

**Nothing was tuned.** No exponent, coefficient or validity bound was adjusted
to reduce these numbers. They are recorded as the measured limit of the
correlation under strong supply-side crossflow with a rounded inlet, and
pinned by `test_rohde_reveals_the_crossflow_limit_of_the_rd_term`.

**Practical bound this puts on use.** For a plenum-fed jet plate -- the #375
case -- `U1/Vi = 0` and none of this applies. For a duct-fed orifice the
correlation is good to roughly `+8%` while `U1/Vi <~ 0.2` (`VHR >~ 25`), and
degrades badly beyond `U1/Vi ~ 0.35` (`VHR <~ 10`), worst for sharp edges.

## Validation targets

Digitisable, in priority order.

0. **M&S Fig. 6, `+` series** -- `(r/d, L/d) = (0, 4.0)`, baseline 0.73.
   **Attempted and not committed.** Automated detection recovered 3 of roughly
   10 markers and picked up annotation text; unlike the open squares, a plus
   encloses no hole, so it cannot be found by the hole detector that made the
   square series reliable, and template matching fires on the frame and on
   text. Needs a manual pass. Committing it would double the validation's
   `t/d` span (0.51 and 4.0).
1. **M&S Fig. 4** -- `C_d` against `U_1/V_i`, compiled data from Rohde, Grimm,
   and Meyfarth and Shine, with the Eq. (17) curve drawn. The scatter here sets
   the harness tolerance for the crossflow term. Visual spread is roughly
   `+/-0.02-0.03` in `C_d`.
2. **M&S Fig. 3** -- `C_d` against `L/d`, Lichtarowitz et al. data with Eq. (13)
   drawn. Visual spread roughly `+/-0.02` in the flat region.
3. **M&S Fig. 2** -- `C_d` against `r/d`, Benedict, Cusick and ASME nozzle data.
4. **Rohde Fig. 8** -- the composite: six faired `t/d` curves (0.51, 1.05,
   1.60, 2.00, 2.83, 4.00), all sharp-edged, at one duct Mach (0.13-0.14),
   against velocity head ratio. The best available test of the `L/d` term
   against data, and it needs the conversion derived above. Surveyed and
   legible; not digitised (the six curves cross, so they cannot be separated
   by y-ordering).
4b. **Rohde Fig. 10** -- **DONE**, digitised by the reviewer. See above.
5. **Rohde Fig. 16** -- single orifice against two in tandem at 1.5 diameters,
   the evidence behind item 33.

Per the README's separation: Figs. 2, 3 and 4 each carry both a **fitted curve**
(already in hand as equations -- no digitising needed) and the **scatter** (what
a tolerance must come from). Only the scatter needs digitising.

---

## Review log

| date | reviewer | outcome |
|---|---|---|
| 2026-09-23 | Claude | **Implemented and validated against the source's own figures.** `orifice::mcgreehan_schotsch` (per-equation stages) + `orifice::Cd_McGreehanSchotsch` + `mcgreehan_schotsch_1988_cd` in Python. Checks A-J moved from the scratch script into 17 gtest cases and 7 pytest cases, each anchored outside the implementation. Two defects were found by the tests, both mine, not the paper's: (i) my Reynolds floor of `Re = 1` let Eq. (8) diverge -- it passes `Cd = 1.0` at `Re = 904` and returned `Cd = 214` at `Re = 1` -- now held at the stated `re_min = 1e4`; (ii) a test asserting the combined `r/d`+`L/d` correction *reduces* the length benefit failed, and the code was right: see item 16a, the paper's prose and its formula point in opposite directions. **The Florschuetz re-check: the Figure 6 scorecard is unchanged at N=242, bias +3.74%, RMS 9.62%, and provably cannot change** -- `C_D` appears zero times in `jet_array_runner.py`, because Fig. 6's `Gc/Gj` is a digitised abscissa. Fig. 6 validates the crossflow bracket *given* `Gc/Gj`; it cannot validate `C_D`. Where `C_D` does act is the predicted `Gc/Gj`: `-3.2%` to `+11.6%` across representative geometries (largest at thin plates, `t/d = 0.5`, at downstream rows), which is `-1.2%` to `+0.4%` in `Nu`. At `t/d = 1` the change is nil by construction, since that is where check H's 0.6% agreement sits. |
| 2026-09-23 | reviewer | **CONFIRMED -- I1 decided: adopt the full chain, crossflow term included, to be validated against the source's own figures.** I2 and I3 follow from that decision (the `U_1/V_i` term is in scope; Eq. (17) stays faithful, bump included). I4 accepted implicitly -- Rohde contributes validity statements and digitisable data, not a correlation. I5 (Wu) deferred, out of scope. Reviewer also asked that a constant-`C_D` path remain available so an implementation or a user can switch back deliberately; see decision D9. |
| 2026-09-23 | Claude | **Scoping correction found while wiring the validation, before any code was written: the crossflow the M&S chain models is not the crossflow Florschuetz models.** See decision D8. `U_1/V_i` is the *inlet/approach* tangential velocity on the SUPPLY side of the plate; Florschuetz's `Gc/Gj` is spent-air crossflow in the impingement channel, on the DISCHARGE side. For Florschuetz's plenum-fed plate the correct evaluation is therefore `U_1/V_i = 0`, which is exactly the configuration check H anchors at 0.6%. The crossflow term is still in scope and still needed -- for a duct-fed orifice, which is Rohde's configuration and a real combaero case -- but it is not what feeds Florschuetz's `beta`. Conflating the two would have applied a discharge-side ratio to a supply-side correlation and biased every `C_D` low. |
| 2026-09-23 | extracted by Claude | UNCONFIRMED -- submitted for review. All three sources read page-image by page-image at 400 dpi. Ten cross-checks run (A-J), all passing; two findings worth the reviewer's attention arrived *through* those checks rather than from reading. **First:** Eq. (17) is non-monotonic in crossflow (item I3). I took this for a transcription error, re-rendered Fig. 4 at 400 dpi to disprove it, and found the bump drawn on the page at the position and height the equation predicts (`0.6344` at `0.085` computed, ~`0.635` at ~`0.1` drawn). The transcription was right and my physical expectation was wrong. **Second:** Eq. (16) is printed `(L/d)' = L/D - r/d` with a capital `D` confirmed by glyph at 400 dpi -- an erratum in the paper, resolved by derivation three ways (item 17a). Also corrected the title recorded in #375 and in the prior session note: "...With Rotation and **Corner Radiusing**", not "...and Crossflow". Rohde is reclassified from a candidate correlation to data plus validity statements -- it contains no equation (item 31) -- but it pays for itself anyway with item 33, which licenses applying a single-hole `C_d` to a jet array above a hole spacing of `1.5d`. No implementation has been written; the gate is I1-I5. |

Change **Status** at the top of this file when reviewed, and record corrections
here rather than silently editing the tables above -- a correction is evidence
about the extraction process, not just about the number.
