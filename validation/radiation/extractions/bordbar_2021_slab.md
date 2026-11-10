# Extraction: Bordbar et al. (2021), 1-D slab radiative source and flux (Figs. 3-5)

**Status: FIDELITY, scored by `validation/radiation/bordbar2021_slab.py`
(#1071).** The data are in `validation/radiation/data/bordbar2021_slab/`:
18 WSGG model curves (`wsgg_curves.csv`, 3609 points) and the 18 sets of
line-by-line squares plotted with them (`lbl_points.csv`, 459 points). The
coordinates are our measurement of the figures; the article is not in git.

## Source, pinned

> Bordbar, H., Coelho, F. R., Fraga, G. C., Franca, F. H. R. and Hostikka, S.
> (2021). Pressure-dependent weighted-sum-of-gray-gases models for
> heterogeneous CO2-H2O mixtures at sub- and super-atmospheric pressure.
> *Int. J. Heat Mass Transfer* 173, 121207.
> doi:10.1016/j.ijheatmasstransfer.2021.121207. CC BY 4.0.
> `sources/radiation/bordbar_et_al_2021_IJHMT_173_121207_pressure_dependent_WSGG.pdf`
> (gitignored). The register entry is in `validation/sources.yaml`.

## What it is

Section 4.2: a 1 m slab between walls at 300 K (taken black here; the
scored residual is the check on that reading), T = 300 K + 1500 K sin^2(pi x/L) (Eq. 18), and three
composition profiles: Case 1 x_H2O = 0.2 sin^2(pi x/L), x_CO2 = x_H2O/2 (Mr
= 2, the rest transparent, Eq. 19); Case 2 x_H2O = 0.1 (1 + sin^2), x_CO2 =
1 - x_H2O (Mr 1/9 - 1/4, Eq. 20); Case 3 x_H2O = (1 + sin^2)/3, x_CO2 = 1 -
x_H2O (Mr 0.5-2, Eq. 21). Fig. 3 is p = 1 atm, Fig. 4 30 atm, Fig. 5 0.5
atm; each has Sr [kW/m^3] and qr [kW/m^2] against x/L for the three cases.
The lines are the new WSGG model, so they test an implementation of it; the
squares are LBL and test the model itself.

## How it was extracted

From the PDF's vector stream (pages 7-8), not from a raster: each panel's
frame, its axis tick marks and its one solid curve path (>= 200 vertices),
and the LBL squares as small closed paths inside the frame. Each axis was
calibrated by a least-squares line through its tick positions against the
printed tick labels. Calibration residual: x at most 4.1e-4 of the unit
axis, y at most 0.37 units (Fig. 5 Case 3 Sr, span 800); five of the 18
panels are exact to 1e-11, the rest within 0.02-0.37 units. Every panel had exactly one solid curve, so no curve was
assigned by colour or judgement. Values are written to 5 significant
digits.

## Scored

Max 1.25 %, mean 0.28 % of each curve's peak |value|, all 3609 points (no
curve above 1.25 %); see the runner's docstring for the quadrature choice.
The LBL squares, through the same solve, give the paper's model errors of
its Table 2 (p = 30 atm Case 2 qr: 10.3 / 3.1 % against the printed 10.0 /
2.82 %).

Not tested by any figure: the 0.1, 2, 5, 10, 20, 50 and 80 atm sets (seven
of ten). Fig. 2 (emissivity charts) is not used: it does not follow Eq. 7
(see the register entry).

## Review log

- 2026-11-09: extracted and scored in scratch (#1071).
- 2026-11-10: committed with the implementation (#1071).
