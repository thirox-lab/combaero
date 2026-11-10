# Radiation gap: what a minimal, sourced combaero radiation capability needs

Design record, 2026-11-02, from screening the radiation collection (validation/sources.yaml, topic radiation). It covers #793, the liner-wall-temperature examples: Lefebvre & Ballal Ch. 8, Topal 2014, the Pretoria T56 work and F101. Today combaero has no radiation at all. `components.py` says "Radiation is not modelled".

## 1. Pieces

A liner wall balance is R1 + C1 = R2 + C2 = K1-2 (Lefebvre Eq. 8.25). combaero already has C1, C2 and K (convection, ThermalWall/WallLayer). Four pieces are missing.

| Piece | Minimal closure | Held source | Validity / caveat |
|---|---|---|---|
| Gas-to-liner radiation R1 | R1 = 0.5 sigma (1 + eps_w) eps_g Tg^1.5 (Tg^2.5 - Tw^2.5) (Eq. 8.7), from Eq. 8.5 with alpha_g/eps_g = (Tg/Tw)^1.5 (Eq. 8.6) | Lefebvre & Ballal 3rd ed. Ch. 8 (held); Andreini 2007 and Topal 2014 reprint it | Primaries: Lefebvre & Herbert 1960 (wanted). eps_w is about 0.7 for Nimonic and 0.8 for stainless steel (book). This is a gray, single-zone, well-mixed estimate. Tg is the mean "effective" gas temperature (Eq. 8.8). |
| Non-luminous eps_g | (A) Reeves: eps_g = 1 - exp[-290 P L (q l_b)^0.5 Tg^-1.5], with P in kPa, q the fuel/air mass ratio and L = 1 (Eqs. 8.9/8.11). (B) Composition-based WSGG: Bordbar 2014 (1 atm) or Bordbar 2021 (0.1-80 atm). | (A) L&B Ch. 8. (B) Bordbar 2014 and 2021 (held, sources/radiation/) | (A) The primary is Reeves 1956, NGTE Memo M285, an industrial gas turbine burning kerosine. The fit range is unknown. It is written in q, so it is implicitly a kerosine-product law and does not carry over to natural gas or H2. (B) 2014: T 500-2400 K, pL 0.01-60 atm m, Mr 0.01-4, 1 atm only. 2021: 300-3000 K, ten discrete pressures, Mr 0.1-4. It needs a p-interpolation choice, and the coefficients are under a NON-COMMERCIAL, no-redistribution licence, so they cannot ship on PyPI without the author's permission. Neither covers Mr -> infinity (H2/NH3 products, Fu et al. 2026). |
| Luminous term | Multiply by the luminosity factor L in Eq. 8.11. L = 336/H^2 (H = fuel H mass %, Eq. 8.15, Lefebvre's engine-data fit) is the default. Alternatives: L = 7.53 (C/H - 5.5)^0.84 (Eq. 8.12, original), 3 (C/H - 5.2)^0.75 (Eq. 8.13), and Kretschmer-Odgers 0.0691 (C/H - 1.82)^2.71 (Eq. 8.14). | L&B Ch. 8 | The data behind L are not held (AFWAL-TR-84-2104, Lefebvre 1984 IJHMT, Kretschmer & Odgers 78-GT-90, all wanted). The four forms disagree, so the shipped form must be named. A soot-volume-fraction route (Liu et al. 2020: kappa_s = 6 pi E(m) f_s / lambda) is out of reach because combaero carries no soot. Geigle 2015 is background only. |
| Mean beam length l_b | l_b = 3.4 V/A (Eq. 8.10, after Fishenden & Saunders 1950). This gives 0.6-0.9 D_L for tubular liners. For annular liners it gives 1.0 D_L (inner) and 1.2 D_L (outer). | L&B Ch. 8; Howell 7th ed. on-line App. B (held) | The textbook value is Le = 0.9 x 4V/A = 3.6 V/A (Howell App. B, Example K.1). The two differ by 6%. Ship 3.4 for fidelity to the #793 examples and document the difference. |
| Liner-to-casing R2 | Gray two-surface enclosure, F = 1: R2 = sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c) (Tw2^4 - T3^4) (Eq. 8.17). Book shortcuts are 0.4 sigma (aluminium casing) and 0.6 sigma (steel casing) (Eqs. 8.18/8.19). | L&B Ch. 8 (Eq. 8.16 is the general form) | Standard enclosure physics; the source only matters for the shortcuts. Use the casing temperature, or T3 as the book does. A_w/A_c is 0.8 for tuboannular liners, otherwise the diameter ratio. |

### Implementation shape (for whoever opens the issue)

The project rules shape the implementation:

- **C++ kernels with (f, J).** R1 is the value; its derivative is dR1/dTw = -1.25 sigma (1 + eps_w) eps_g Tg^1.5 Tw^1.5. Tg enters through eps_g.
- **Derivatives with respect to Tg.** d eps_g/dTg is analytic for Reeves.
- **The wall's two radiation terms.** R2 is cold-side and R1 hot-side, attached to ThermalWall like the existing convective sides.
- **Inputs.** Pressure, Tg, q (or composition), l_b from the passage geometry, L from the fuel's H mass fraction, and eps_w and eps_c as user inputs with the book's defaults.

## 2. What can falsify it (VALIDATION_POLICY split)

1. **Fidelity: Lefebvre & Ballal Ch. 8 worked examples (held).** These are printed-quantity fidelity checks. The uncooled liner has Tw1/Tw2 = 1640/1603 K, the film-cooled case 1283/1265 K and the TBC case 1304/1236/1220 K. Re-solving the printed chain already reproduces them within 1 K (register entry). The test is that a combaero network with the new R1/R2 kernels gives the same values, using the book's Tg = 2280 K as an input. This exercises Eqs. 8.7, 8.9/8.11, 8.10 and 8.17-8.19 together. It cannot test the physics of eps_g, only our implementation.
2. **Kernel cross-check, accuracy only.** At 1 atm, compare Reeves eps_g with Bordbar 2014 on kerosine products for several q and l_b. Two further one-number checks from Howell (Alberti spreadsheet): Example K.1 (eps_g 0.1451: ethane products, 1853 K, Le 1.44 m) and Example F.1 (0.401, Le 2.4 m). Both models are model-to-model, so any gap is labelled as a model difference, not as error.
3. **Cross-source accuracy against measured radiant flux: Humenik, Claus & Neely (NASA TM-83436 and TM-83538, wanted, public).** The rig is a tubular-can combustor tested at 0.34-2.07 MPa, 533-700 K and F/A 0.008-0.021. Two fuels were run: Jet A (13.9% H) and ERBS (12.9% H). Total and spectral radiometers sit at three stations, and liner temperatures were measured. This is the first case that tests L against fuel H content and pressure on a different rig. Rudey & Graham 1984 (held) show it reduced to radiant flux per station (idle 269 kPa, cruise 700 kPa, take-off 1740 kPa). Usable only if the TMs print the can geometry and zone gas temperatures; unchecked.
4. **F101 (Lefebvre Fig. 8.5), via AFWAL-TR-84-2104 (held since 2026-11-03).** F101 is one of the six USAF programs Eq. 8.15 was fitted over, not the fit itself; the fit is in AFWAL-TR-83-2004 (wanted). AFWAL-TR-84-2104 compares the measured F101 wall temperature only with uncooled calculations and prints no geometry or per-point conditions, so F101 is neither a fidelity nor an accuracy case yet. Gleason et al. 1979, AFAPL-TR-79-2018 (wanted), should hold what is missing.
5. **Weak cases.** Ballal & Lefebvre SME4 Figs. 5/6: thermal-paint ranges, gas T from an assumed efficiency. Topal 2014: measured paint and thermocouple temperatures, but the geometry, FAR and L are not printed. Treat both as qualitative.

## 3. Missing: wanted-list proposals

Add these:

- **Humenik, F. M., Claus, R. W., Neely, G. M. (1983). Parametric study of flame radiation characteristics of a tubular-can combustor. NASA TM-83436.** Public on NTRS: https://ntrs.nasa.gov/api/citations/19830020966/downloads/19830020966.pdf. It settles validation case 3: measured radiant flux and liner temperature against P, T3, F/A and fuel H.
- **Claus, R. W., Neely, G. M., Humenik, F. M. (1984). Flame radiation and liner heat transfer in a tubular-can combustor. NASA TM-83538 / AIAA 84-0443.** Public on NTRS: https://ntrs.nasa.gov/api/citations/19840005120/downloads/19840005120.pdf. It adds liner heat loads and flame emissivity derived from spectra, which is the direct check of eps_g.
- **Alberti, M., Weber, R., Mancini, M. (2018). Gray gas emissivities for H2O-CO2-CO-N2 mixtures. JQSRT 219, 274-291 (doi:10.1016/j.jqsrt.2018.08.008).** Add with its companion: **(2020) New formulae for gray gas absorptivities of H2O, CO2 and CO. JQSRT 255, 107227.** These give a pressure-capable (0.1-100 atm, 300-3000 K) total emissivity and absorptivity in closed form, with a public spreadsheet, and they cover single gases (pure H2O, so the H2/NH3 limit). They are the likely shippable alternative to Bordbar 2021's non-commercial coefficients; check the licence of the fit coefficients. Elsevier, possibly open access.
- **Kretschmer, D., Odgers, J. (1978). A simple method for the prediction of wall temperatures in a gas turbine combustor. ASME 78-GT-90.** Source of Eq. 8.14, and a second liner method with data. ASME Digital Collection.
- **Stuttaford, P. J., Rubini, P. A. (1998). Assessment of a radiative heat transfer model for gas turbine combustor preliminary design. J. Propulsion and Power 14(1), 66-73.** This is prior art for exactly what combaero would do (radiation inside a combustor network solver), and it has an annular-liner temperature comparison. AIAA (paywalled). Found as a citation in the deleted Malalasekera review.
- Optional, priority later:
  - Taylor, P. B., Foster, P. J. (1974). The total emissivities of luminous and non-luminous flames. IJHMT 17, 1591-1605. A soot-plus-gas WSGG for oil flames, an alternative luminous closure. Elsevier.
  - Leckner, B. (1972). Spectral and total emissivity of water vapor and carbon dioxide. Combust. Flame 19, 33-48. The classic closed-form total emissivity with a pressure correction and no licence encumbrance. Elsevier.
  - Reeves, D. (1956). NGTE Memo M285. The primary of Eq. 8.9 and its fit range. UK National Archives (DSIR/AVIA series), likely unobtainable.

Keep these, since none was obtained in this batch:

- Lefebvre & Herbert 1960 (Proc. IMechE 174).
- Lefebvre 1984 (IJHMT 27, 1493-1510). Liu et al. 2020 and Howell App. K both cite it as the gas-turbine flame-radiation reference.
- AFWAL-TR-84-2104 (DTIC ADA151464, browser download).

Remove: nothing.

## 4. NASA tubular-can case (screened)

Both TMs were screened from page images on 2026-11-02. Draft entries are in `entries.yaml`. The conclusion is that validation case 3 in section 2 is weaker than hoped.

- **No geometry.** Neither report prints the liner diameter, length, louvers/slots, holes, annulus, air split or pressure drop. The only dimensions are the 25.4 cm housing and the instrument stations: total gauges at x = 12.7, 18.5 and 28.7 cm, spectral ports at 7.4, 22.1 and 34.0 cm. The thermocouple layout implies 10 louvers, undimensioned. The liner carries a TBC on Hastelloy, of unstated thickness. **It cannot be a combaero liner network case.**
- **Conditions are complete.** P3 0.34-2.07 MPa, T3 533/616/700 K, the airflow per point (constant flow parameter, U_ref 15 m/s on the housing) and f/a 0.016 (sweep 0.0087-0.0212) are all printed. So are Jet A (13.9 % H, 17.2 % aromatics) and ERBS (12.9 % H, 28.8 %).
- **Data are plotted only, without uncertainty.** That covers the total radiant flux per zone, the outer-wall (T_liner - T3) per zone, and 12 flame-emissivity points: Fig. 7 of TM-83538, Jet A, 533 K, f/a 0.016. There is no gas temperature, no zone f/a and no l_b.
- **Usable as an accuracy test of L, in ratio form.** ERBS/Jet A flux at equal P3, T3, f/a and nozzle, where the TM says Tg is unchanged between fuels. Measured: +56 % (primary, 0.34 MPa) falling to +5 % (2.0 MPa); none in the tertiary zone. A screening estimate (not a record) with Fig. 7's emissivities gives:

  | L form | Predicted ERBS/Jet A ratio (0.34 / 2.0 MPa) |
  |---|---|
  | L = 336/H^2 | 1.14 / 1.06 |
  | The C/H forms | 1.31-1.50 / 1.13-1.18 |

  The 336/H^2 form under-predicts the low-pressure fuel effect. The C/H forms come closer at low pressure but over-predict at 2 MPa. No form fits both ends. This is the cross-source (same lab, different rig to F101) accuracy evidence for picking the L form; it is not fidelity.
- **Absolute flux or eps_g against P3.** This is possible only with an assumed l_b (liner diameter) and Tg, so it carries assumption labels.
- **The TM's own model.** TM-83538's liner model is NREC 1980 plus a Juhasz-Marek film with C_m 0.05, and alpha_w 0.4 for the TBC. Its non-luminous formula is misprinted, so it is a method check only.
- **Next leads (wanted).** These may hold the geometry:
  - Srivatsa 1982, NASA CR-167930, which predicted this rig (NTRS 19820021901; scripted downloads truncate, fetch in a browser).
  - Claus 1981, NASA TP-1722.
  - NREC 1344-2 (1980).

## 5. Decision the user has to make

Which non-luminous kernel ships?

- **Reeves (Eq. 8.9) alone.** Simplest. It matches the #793 examples, but it is kerosine-only and has no stated range.
- **Reeves plus a composition-based option.** The option would be Bordbar 2014 (shippable? its licence is unchecked; the 2014 paper prints Table 1 in the article) or Alberti (once obtained). Bordbar 2021 cannot ship without the author's permission.

The luminous factor L can only be the Lefebvre family until a soot model exists.

## 6. Outcome: Bordbar 2021 shipped as the opt-in option (#1071, 2026-11-10)

The licence note in sections 1 and 5 was wrong: the article is CC BY 4.0 and Appendix A (supplement mmc1) holds every constant; only the mmc2 MATLAB/.dat files carry the non-commercial notice, and they are not used. The user chose Eq. 8.11 as the default and Bordbar 2021 opt-in (`GasRadiation(emissivity_model="bordbar2021")`).

- **Formulation.** Eqs. 7-11 as printed; constants from Appendix A only (`include/wsgg_bordbar2021_data.h`). Between the ten sets eps is interpolated linearly in ln P, each set evaluated at the actual pressure path. Rejected: interpolating the coefficients (gray gas i is a different band in each set; the 0.5 atm set's gas 4 is its weakest) and evaluating each set at its own tabulated pressure (neighbouring sets then differ 40-75% at one state, against 1-4% at the actual path (DLR points), so the result leans on the interpolation rule). Linear-in-P and ln-P weights agree to 1e-4 at the DLR points; ln P was taken because the sets are spaced roughly geometrically.
- **Invented, not forced by the paper:** the interpolation rule; clamping T to 300-3000 K, Mr to 0.1-4 (no CO2 is Mr = 4, an extrapolation for H2/NH3 products) and eps to [0, 1], with zero partials where clamped; the end set outside 0.1-80 atm.
- **Fidelity.** Figs. 3-5 slab curves (1, 30, 0.5 atm; 3609 points from the PDF's vector paths) re-solved with the shipped kernel: max 1.25%, mean 0.28% at 16 Gauss-Legendre ordinates per hemisphere, converged (1.42 / 0.44% at the paper's 6, where Gauss-Legendre is not yet converged). The bound catches gross errors and a 5% change of a dominant gas's weight, not a few-% change of a weak gas; the transcription itself was checked against the supplement's own files (identical to 2.5e-16). Seven sets are untested by any figure.
- **Jacobian.** The emissivity enters the ThermalWall residuals, so its T, P, composition and path partials are analytic (C++) and checked against central differences at kernel, coupling and network-relay level; each was falsified (dropping the ln P weight's P term, dk_dY or dk_dT turns the relay or coupling tests red). Radiation off and the Eq. 8.11 path are bit-identical to before (hash of the effusion runner and ThermalWall network solutions).
- **Accuracy (cross-source).** Behrendt 2008 effusion rig: eps_g 0.132 / 0.198 at 5 / 10 bar against Eq. 8.11's 0.137 / 0.262. The score does not improve (flat plate +9.2% per vertex against +8.9%), and the 5 -> 10 bar trend closes less (0.58x against 0.84x); see `validation/cooling/extractions/behrendt_2008_effusion_eta_tot.md`.
- **Tests:** `python/tests/test_wsgg_bordbar2021.py`, `tests/test_wsgg_bordbar2021.cpp`, the `bordbar2021` cases of `python/tests/test_liner_radiation_793.py`, `PINNED_RAD_BORDBAR` in `python/tests/test_effusion_network_runner.py`.
