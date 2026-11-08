# Extraction: Behrendt, Lengyel, Hassa and Gerendas (2008), effusion eta_tot in a reacting 5-10 bar flow

**Status: SCORED (accuracy, cross-source). Extracted 2026-11-08 from the PDF
vector stream; the user confirmed it on the overlay review ("looks
perfect").** The data are in `validation/cooling/data/behrendt2008/`: 9
distinct `eta_tot(x/d)` curves (Figs. 6, 9, 10) and the 3 Fig. 8 heat-flux
series. The 9 eta curves are scored by `validation/cooling/effusion_network_runner.py`
through `EffusionPlateElement` (set `andrews_1986_effusion_internal`, the
element's coolant side), at both gas-side baselines; see "Scored" below. The
Fig. 8 series stay `scores: null`.

## Source, pinned

> Behrendt, T., Lengyel, T., Hassa, C. and Gerendas, M. (2008).
> Characterization of advanced combustor cooling concepts under realistic
> operating conditions. *ASME Turbo Expo 2008*, GT2008-51191, Vol. 4,
> 1801-1814. doi:10.1115/GT2008-51191.
> `sources/heat_transfer/film/Behrendt_Lengyel_Hassa_Gerendas_2008_GT2008-51191_effusion_reacting_high_pressure.pdf`
> (gitignored, copyrighted). The register entry is in `validation/sources.yaml`.

Rig paper: Behrendt and Hassa (2006), ICAS 2006 paper 152 (held; the 2008
IMechE journal paper has the same content). The authors' wall model, which
gives their M, I and Table 4: Martiny, Schulz and Wittig (1997), 97-GT-329
(held). Page numbers below are Vol. 4 pages: PDF page n is page 1800 + n.

## What it is

The hot gas comes from a kerosene combustor. It flows past a 100 x 80 mm
effusion sample in the side wall of a rectangular duct. The rig runs at
p_h 5 bar (TC1) or 10 bar (TC2), with near-wall gas at 1567 K and 1532 K.
The wall temperature is measured by IR through the opposite window. The
sample is Nimonic 75 and conducts, so the quantity is the **total (overall)**
effectiveness, not the adiabatic film effectiveness. The authors say so
(p. 1807). Eq. (2):

    eta_tot(x/d, z/d) = (T_h(x/d=-10, z/d) - T_w(x/d, z/d)) / (T_h(x/d=-10, z/d) - T_c)

- **T_h** is the near-wall gas temperature at x/d -10, y/d 4, measured by
  thermocouple (ICAS 2006). It is not a bulk temperature.
- **T_c** is the nominal coolant temperature, measured y/d ~15 behind the
  sample.
- **Lateral average:** the profiles are averaged over -17.1 < z/d < 8.6,
  which is three lateral pitches.

In `conventions.py` the source is declared as `eta_overall`, with a comment
on the T_h reference.

**Uncertainty.** The authors give d eta_r,tot <= 3% (p. 1807). The subscript
r means RELATIVE: their "2...3%" dp effect in the same notation is
0.790 -> 0.805 on Fig. 9. This 3% is a worst-case bound built from three
terms: T_h fluctuation of +/-25 K, the thermocouple junction depth (7 K) and
emissivity scatter (+/-12 K). It is not one sigma. The schema's band is
relative too, so every eta series carries `uncertainty: 0.03` unchanged. The
number is a measurement statement, constant over the source, and was never
derived from a model's error (#389).

## Inputs table

| item | value | where |
|---|---|---|
| samples | E20 (beta 20 deg), E30 (beta 30 deg); no compound angle | Table 1, p. 1804 |
| hole d | **~0.70 mm, inferred**. l/d 143 with l = 100 mm gives 0.699; b/d 114 with b = 80 mm gives 0.702 | p. 1803 + ICAS 2006 |
| d check | Hausen Eq. (5) at d = 0.70 mm against Table 4's h_e: Nu 16.3 (table) vs 16.8 (Eq. 5) at Re_e 4600; 34 vs 35.3 at 9650 | p. 1810-1811 |
| wall | t/d 2, NIMONIC 75 (2.4951); L/d = (t/d)/sin(beta): 5.85 (E20), 4.0 (E30), derived | p. 1804 |
| pattern | sx/d 5, sz/d 8.6, rows staggered, repeats every 5 rows | p. 1804 |
| upstream | no starter film; BL bleed at dx/d -21 | pp. 1803-1804 |
| hot gas | TC1 5 bar, T_h 1567 K, u_h 29.7 m/s; TC2 10 bar, 1532 K, 29.1 m/s (near-wall, laterally averaged) | Table 2 |
| turbulence | Tu 19-24% | ICAS 2006 |
| coolant | T_c 470-580 K, dp_r,c 2-4%; M and I per legend | Table 3, legends |
| duct | A = 80e3 mm2; combustor ~300 mm upstream of the sample | ICAS 2006, p. 1803 |

**Caveats on the inputs.**

- **M and I are not measured.** They come from the authors' 1-D wall model
  (Martiny et al. 1997), with cD from Martiny's 1999 dissertation. The model
  puts the hole exit 40-70 K above T_c.
- **About 40% of the preheated coolant cooled the module, not the sample**
  (footnote 3, p. 1811). The split was measured only under cold conditions,
  so the sample's mass flow, and hence M, is approximate.
- **T_h is not uniform laterally:** 1436-1570 K across the duct (Fig. 3).
- **Test-case labels differ between papers.** This paper's TC1 and TC2 are
  not the ICAS paper's TC1 to TC4.

## Figures extracted

| Fig. | page | content | committed |
|---|---|---|---|
| 5 | 1809 | E20 vs E30, TC1, dp 3%, Tc 540 K | none of its own; both curves are sweep members (below) |
| 6 | 1810 | E30, TC1, Tc 470 / 540 / 580 K | 3 curves, to x/d 58 |
| 8 | 1811 | hot-side q vs dp_r,c: TC1 x/d 26 and 74, TC2 x/d 26 | 3 x 3 markers, **biased, not for scoring** |
| 9 | 1812 | E20, TC1, Tc 540 K, dp 2 / 3 / 4% | 3 curves, to x/d 50 |
| 10 | 1812 | E30, 5 and 10 bar, dp ~2 / 3% | 3 curves, to x/d 58 (5 bar dp 3% is the Fig. 6 Tc 470 K copy) |

Figs. 4 and 7 are contour maps and were not extracted. The I20/I30
impingement-backed samples appear only in Fig. 4.

## Method and calibration evidence

1. **Vector paths, not pixels.** PyMuPDF `get_drawings()` on each page, in
   CropBox coordinates. Each eta curve is one stroked polyline with no gaps
   and x-monotonic vertices. The vertices are the plotted points: spacing
   0.340 x/d on E30 and 0.367 x/d on E20, with a few 2x steps where a
   collinear vertex was dropped. The Fig. 8 data are the centres of the
   4.52 pt square markers. The connecting polyline's vertices coincide with
   them to < 0.001 pt.
2. **Calibration.** Each axis is a linear least-squares map from pt to data
   units, fitted on every labelled major tick (labels from the page text
   layer):

   | Fig. | x ticks | max x residual | y ticks | max y residual | tick-spacing max dev (x / y) |
   |---|---|---|---|---|---|
   | 5 | 6 | 1.0e-3 x/d | 4 | 3.2e-7 | 0.0097 / 0.0077 pt |
   | 6 | 6 | 8.7e-15 x/d | 4 | 3.2e-7 | 0.0008 / 0.0009 pt |
   | 8 | 7 | 7.3e-6 % | 7 | 1.5e-4 W/cm2 | 0.0007 / 0.0010 pt |
   | 9 | 6 | 7.5e-5 x/d | 7 | 3.6e-7 | 0.0009 / 0.0008 pt |
   | 10 | 6 | 2.0e-4 x/d | 4 | 3.6e-7 | 0.0008 / 0.0005 pt |

   The tick-spacing check uses all 16-31 ticks per axis, minor ones
   included. Uniform spacing confirms linear axes. It also rules out a tick
   label assigned to the wrong tick.
3. **Overlay.** The page was rasterised independently by poppler
   (`pdftoppm -cropbox`), and the extracted points were drawn over it for
   the user's review.
4. **CropBox trap.** PyMuPDF reports drawings in CropBox coordinates. Poppler
   renders the MediaBox unless told `-cropbox`. On these pages the CropBox
   is the MediaBox (595 x 842 pt) trimmed by 24.98 pt top and bottom, so an
   overlay rendered without `-cropbox` is offset by 24.98 pt vertically and
   looks like a calibration error. It is not one. The fix is `-cropbox`, not a
   shift of the data.
5. **Clipped last vertex dropped.** Every eta polyline ends exactly on the
   right frame (x/d 50.0006 / 58.0001). That vertex is the plotting
   program's clip point, not a station, so it is not in the CSVs and is
   excluded from the plateaus.
6. **Fig. 8 abscissae.** The markers land on the run settings 2.0/3.0/4.0
   (2.1 at TC2) to 5e-6. This checks the x calibration independently.

## The duplicate cross-check

Three curves are printed twice, in different figures with independent
calibrations. The two copies of each were compared on the first copy's
abscissae (linear interpolation; the last 0.5 x/d is excluded):

| copy A | copy B | n | max \|d eta\| | committed |
|---|---|---|---|---|
| Fig. 5 E20 | Fig. 9 dp 3% | 85 | 1.05e-5 | Fig. 9 |
| Fig. 5 E30 | Fig. 6 Tc 540 K | 93 | 1.0e-5 | Fig. 6 (drawn to x/d 58; Fig. 5 stops at 50) |
| Fig. 6 Tc 470 K | Fig. 10 5 bar dp 3% | 116 | 9e-7 | Fig. 6 |

Agreement to 1e-5 tests three things together: the series assignment, the
colour mapping and both y calibrations. Only one copy of each curve is
committed, and the rule is the earliest one-parameter sweep figure (6, 9,
10). Fig. 5 only compares two sweep members, so it has no series of its own
and no redraw. The Fig. 10 redraw lacks its 5 bar dp 3% curve.

## Plateau table (x/d 18.1-49.8)

The authors exclude x/d < 16, which an engine's starter film would dominate,
and the curves start at x/d 18.1. Each mean is the trapezoid average over
the vertices up to x/d 50, without the clip vertex. The window is 18.1-49.8
for every curve, including those drawn to 58.

| Fig. | sample | p_h [bar] | dp_r,c [%] | T_c [K] | M | I | n | mean | min | max |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | E30 | 5 | 3.0 | 470 | 5.6 | 10.5 | 94 | 0.7780 | 0.7697 | 0.7882 |
| 6 | E30 | 5 | 3.1 | 540 | 5.3 | 10.8 | 94 | 0.7266 | 0.7168 | 0.7391 |
| 6 | E30 | 5 | 3.0 | 580 | 5.0 | 10.3 | 93 | 0.6883 | 0.6798 | 0.7008 |
| 9 | E20 | 5 | 2 | 540 | 4.1 | 6.5 | 85 | 0.7900 | 0.7824 | 0.8053 |
| 9 | E20 | 5 | 3 | 540 | 5.2 | 10.3 | 86 | 0.7996 | 0.7912 | 0.8152 |
| 9 | E20 | 5 | 4 | 540 | 6.1 | 14.2 | 87 | 0.8051 | 0.7957 | 0.8222 |
| 10 | E30 | 5 | 2.0 | 470 | 4.5 | 6.7 | 94 | 0.7681 | 0.7610 | 0.7782 |
| 10 | E30 | 10 | 3.0 | 480 | 5.8 | 11.7 | 94 | 0.7529 | 0.7392 | 0.7674 |
| 10 | E30 | 10 | 2.1 | 480 | 4.8 | 8.0 | 94 | 0.7358 | 0.7212 | 0.7510 |

Fig. 8 (biased): at x/d 26, q rises from 45.5 to 53.7 to 58.1 W/cm2 at TC1
and from 65.6 to 68.7 to 72.9 W/cm2 at TC2, over dp 2 (2.1) / 3 / 4%. At
x/d 74 (TC1) it stays flat at 26.0-28.0 W/cm2.

## What it can and cannot settle

**It can settle:**

- **Accuracy, not fidelity.** It is a cross-source ACCURACY set for
  `EffusionPlateElement`'s `eta_overall` at reacting 5-10 bar conditions.
  It must be labelled cross-source in any report. The authors offer no
  correlation, so there is nothing whose fidelity this could judge.
- **Four one-parameter trends.** Each is a plateau mean from the table above:

  | trend | condition | eta_tot |
  |---|---|---|
  | hole angle 20 vs 30 deg | TC1, dp 3%, Tc 540 K | 0.800 vs 0.727 |
  | coolant temperature 470 / 540 / 580 K | E30, TC1, dp 3% | 0.778 / 0.727 / 0.688 |
  | pressure drop 2 / 3 / 4% | E20, TC1, Tc 540 K | 0.790 / 0.800 / 0.805 |
  | p_h 5 -> 10 bar | E30, dp 3% | 0.778 -> 0.753 |
  | p_h 5 -> 10 bar | E30, dp ~2% | 0.768 -> 0.736 |

  The two p_h rows are at Tc 470 K and 480 K respectively. The authors put
  the effect of that 10 K difference at about 1.3% relative (p. 1811).

  A model that gets the level right and a trend's sign wrong is still wrong.
- **Never a tuner** (docs/VALIDATION_POLICY.md). Matching this rig is the
  user's job through `Nu_multiplier`. The harness will not score a
  multiplier fitted here, recommend one, or apply one.

**It cannot settle:**

- **gas_augmentation (#476).** The data contain no adiabatic wall
  temperature, so the gas-side h and eta_aw cannot be separated. Like
  Jackowski 2021, the set is a check on the coupled wall, not on its
  gas-side closure. Table 4's h_h values are the authors' model output
  (Martiny superposition), not a measurement.
- **Fig. 8 is not for scoring.** The sensor replaces two cooling holes and
  biases the heat exchange toward the coolant side (p. 1811).
- **The hole Cd.** The coolant flow split is approximate (footnote 3).

**A regime the element does not model (follow-up, not this PR).** Martiny
et al. (1997) define the enthalpy-to-conduction ratio K3 = rho u cp d / k_w.
It is about 6 on these plates. That is below the ~10 where coolant
saturation inside the hole starts to matter (Martiny Fig. 14): the coolant
heats toward the wall temperature along the hole, which caps the heat the
hole can remove. `EffusionPlateElement` has no coolant capacity limit and
takes the sink at the supply temperature, so on this set it may over-predict
the internal cooling. Test that before reading any level error as a
gas-side gap.

## Scored: EffusionPlateElement, both gas-side baselines (accuracy, cross-source)

`effusion_network_runner.py` builds the rig as a network (gas MassFlowBoundary
-> MomentumChamberNode, A = 80e3 mm2 -> exit; coolant PressureBoundary at
p_h (1 + dp_rc), T_c -> the plate -> chamber) and reads the plate's own
`eta_overall`. The panel is the scored window x/d 18-50 (32 d x 80 mm): the
element is homogenised and predicts one eta per panel, so it is placed where
the measurement is; vertices beyond x/d 50 are another window and come back
unscored. `gas_augmentation`, `internal_Nu_multiplier` and the chamber's
`Nu_multiplier` are 1.0, `gas_film` none. Nothing is fitted. Pinned by
`python/tests/test_effusion_network_runner.py`.

- **chamber** (the default): the chamber's SmoothModel(gnielinski), fully
  developed, Dh = sqrt(4A/pi) = 319 mm.
- **flat_plate**: `gas_baseline="flat_plate"`, the boundary layer restarted at
  the bleed (dx/d -21), so x in [39 d, 71 d]; CR-2648's St_o averaged over it.
  Re_x 1.6-2.9e4 (5 bar) and 3.2-5.9e4 (10 bar) are below the Stanford data's
  1.5e5, so these rows are out-of-domain by the kernel's own flag.

Plateau means (x/d 18-50), predicted at gas_augmentation 1.0:

| curve | measured | chamber | err | flat plate | err | h0 chamber / flat [W/m2K] |
|---|---|---|---|---|---|---|
| E20 dp 2% Tc 540 K | 0.790 | 0.926 | +17.2% | 0.860 | +8.8% | 94 / 193 |
| E20 dp 3% Tc 540 K | 0.800 | 0.934 | +16.8% | 0.874 | +9.3% | 94 / 193 |
| E20 dp 4% Tc 540 K | 0.805 | 0.939 | +16.7% | 0.883 | +9.7% | 94 / 193 |
| E30 dp 3% Tc 470 K | 0.778 | 0.924 | +18.8% | 0.856 | +10.0% | 94 / 193 |
| E30 dp 3.1% Tc 540 K | 0.727 | 0.925 | +27.3% | 0.857 | +17.9% | 94 / 193 |
| E30 dp 3% Tc 580 K | 0.688 | 0.924 | +34.2% | 0.856 | +24.3% | 94 / 193 |
| E30 dp 2% Tc 470 K | 0.768 | 0.915 | +19.1% | 0.840 | +9.4% | 94 / 193 |
| E30 10 bar dp 3% Tc 480 K | 0.753 | 0.914 | +21.4% | 0.839 | +11.4% | 163 / 334 |
| E30 10 bar dp 2.1% Tc 480 K | 0.736 | 0.906 | +23.1% | 0.824 | +12.0% | 163 / 334 |

Mean +21.6% (chamber) and +12.5% (flat plate); per vertex (821 points, the
scorecard rows) +21.7% and +12.6%. No curve is inside the authors' 3%
band on either baseline. The restarted boundary layer halves the miss and does
not close it.

Trends, plateau differences (measured / chamber / flat plate):

| trend | measured | chamber | flat plate |
|---|---|---|---|
| angle 20 -> 30 deg (TC1 dp 3% Tc 540 K) | -0.073 | -0.010 | -0.017 |
| T_c 470 -> 580 K (E30 TC1 dp 3%) | -0.090 | +0.000 | +0.000 |
| dp 2 -> 4% (E20 TC1 Tc 540 K) | +0.015 | +0.013 | +0.023 |
| p_h 5 -> 10 bar (E30 dp 3%) | -0.025 | -0.010 | -0.017 |
| p_h 5 -> 10 bar (E30 dp ~2%) | -0.032 | -0.009 | -0.016 |

Every sign is right except T_c, which both miss entirely: the two-resistance
closure is T_c-invariant to first order. The flat plate raises the gas-side
weight, so the internal-side trends (dp) grow past the measurement (1.5x)
and the pressure trend closes from 0.3-0.4x to 0.5-0.7x. The angle effect
(the authors: 46% less wall-normal jet momentum at 20 deg) stays at 0.1-0.2x:
the element has no angle-dependent gas side.

Sensitivity (reported, never chosen; `sensitivity()`):

- **Start length moves the score.** Bleed at 10.5 d / 21 d / 42 d: eta 0.850 /
  0.856 / 0.864 (E30 Tc 470 K); from the sample's own leading edge 0.843.
- **Duct area (80e3 vs 8e3 mm2, Dh 319 vs 89 mm) moves only the chamber
  baseline** (-0.02 in eta); the flat plate does not see Dh (identical to 1e-9).
- **Whole sample as one panel** (x/d 0-143): flat-plate eta +0.009 against the
  window, a thicker average boundary layer.

**What remains (not fixed here).** A diagnostic inversion of the element's own
wall (never applied) puts the gas-side h the plateaus need at 1.6-1.8x this
flat plate on seven curves and 2.3 / 2.7x on the two hot-coolant E30 curves.
That is the job `gas_augmentation` (the user's) and the missing near-wall
mixing physics share: a T_c- and angle-dependent film/mixing temperature at
M 4-6, I 6-14, Tu ~20%. Two further gaps sit on the coolant side and are
open: Andrews' throat Nu runs ~1.8x the authors' in-hole Eq. (5), and the
element's coolant flow ~1.4-1.5x the printed M (the printed M are themselves
approximate, footnote 3). Neither is a tuning target.

## Review log

- 2026-11-08: extraction, overlay review page; user: "looks perfect".
  Promoted with one copy per duplicated curve and without the clip vertex.
- 2026-11-08: scored through EffusionPlateElement (effusion_network_runner),
  chamber and flat-plate gas baselines; the user chose the flat plate as the
  baseline to add ("flat plate approach now, it is what we hold").
