# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.11.1] - 2026-11-07

### Fixed

- **hybr no longer turns a nearly singular Jacobian into NaN iterates (#1027).** This is the hybr twin of #1020. MINPACK's hybrj dogleg treats a pivot as singular only when it is exactly zero.
  - **Before:** on stress effusion 105, the Jacobian hybr requested at its 48th evaluation had a pivot of 7.6e-309. An effusion plate's crossflow Cd had underflowed, and its row was the only one that set chamber ln__b1's pressure level. The Gauss-Newton step overflowed, the dogleg step and the trust radius became NaN, and hybr's last 8 evaluations were all NaN. The LM fallback recovered, so the result was right, but hybr's tail was wasted. Effusion 103 had the same 8 NaN iterates.
  - **Now:** hybr gets the same Jacobian LM has had since #1020, with entries below machine epsilon times their column's largest set to zero. The dogleg's zero-pivot guard then applies. Stress effusion 0:200 has no NaN iterates (16 before), and every case still converges. The model's residuals and Jacobian are unchanged.

- **The zero-fuel smoothing no longer under-burns at every equivalence ratio (#1037).** `smooth_phi0` (on by default in the Python `cb.complete_combustion`, the network combustor and the fuel controller's T_ad) smooths the O2 demand a at zero fuel.
  - **Before:** the softplus `(ln(1 + e^(k0 a)) - ln2) / k0` tends to `a - ln2/k0`, so 3.5e-5 mol O2 per mol of mixture stayed unburnt at every lean phi: T_ad 0.33-0.58 K low (CH4/air at 300 K: -0.47 K at Phi 0.02, -0.34 K at Phi 1; H2 -0.48 to -0.58 K). The Python and C++ README quick starts printed 2326.66 K and 2326.99 K for the same mixture. The products Jacobian at exactly zero fuel also missed the burnt share (one-sided slope 100% off).
  - **Now:** `n_req = (ln(1 + e^(k0 a)) - ln2 sech(k0 a)) / k0`. It keeps n_req(0) = 0, the same slope 1/2 at zero fuel and the same transition width (about 1/k0), is infinitely smooth for a > 0, and equals a to round-off once k0 a > ~40. T_ad with the smoothing on is now identical to the unsmoothed value from Phi 0.02 to 2 (CH4, H2, n-heptane, NH3). k0 and every signature are unchanged.

- **The GUI no longer turns an infinite number into a null that the solver refuses (#1038).** A network file with an overflowing literal such as `"pressure_ratio": 1e999` parses as Infinity. Typing `Infinity` or `1e999` into a numeric field also gave Infinity, because the field checked the value while typing but not when it lost focus. Save then wrote `null`, the reloaded file passed the GUI's own checks, which read null as the default, and Solve refused it with "Input should be a valid number".
  - **Load** refuses a file that holds a non-finite number, names the node and field in the note above the canvas, and leaves the canvas as it was.
  - **Numeric fields** reject a non-finite value when they lose focus and go back to the previous value, as they do for text.
  - **Save** refuses to write a non-finite number and says which one, so it never writes `null` for a number.
  - **The pre-Solve check** refuses an explicit `null` in a field the backend types as a plain number with a default: compressor/turbine pressure ratio and efficiency, the tee's angle, loss factor, design share and psi, the ejector's recovery efficiency, and the discrete loss's zeta and zeta0. Previously it read null as the default. A python test keeps this list in step with the backend schema.

- **The Inspector shows an in-pipe orifice's refusal before Solve, whichever pipe causes it (#1039).** An orifice whose bore is wider than its downstream pipe (the one that sets beta in reverse flow, #788) was refused at Solve, but the Inspector compared the bore only with the upstream pipe and showed no note. It now shows the solver's own reason from the resolve preview whenever the orifice is unresolved, naming the pipe. That reason replaces the old upstream-only note, so the upstream case is reported once, not twice.

## [0.11.0] - 2026-11-07

**Upgrading from 0.10.** No public function or argument was removed or renamed. 0.11 adds `FannoChannelFlow::M_exit_root`, `OrificeElement.closed` and `MomentumChamberNode.pipe_area_warnings()`, and it is stricter about geometry:

- **Now refused**, each with a message naming the input:
  - a bore or hole area wider than its pipe or chamber;
  - an explicit area-change, tee or ejector area that contradicts its pipe;
  - overlapping holes;
  - a zero or negative duct or node area;
  - a Dh larger than the area allows;
  - an infinite pressure ratio;
  - an unknown combustor `method`;
  - an energy fraction at or below -1;
  - a `ThermalWall` on a node with no wetted surface.
- **Now accepted:**
  - an area change with F0 == F1, as a pass-through;
  - a zero-area restriction, as a closed valve.
- **Results that change:**
  - a node with no surface no longer takes its flow area as one, so it reports no htc (NaN, "—" in the GUI);
  - Mach below the 1 Pa drive floor is now the flow's;
  - a chamber narrower than its pipe now warns.
- **sympy** is no longer installed with `combaero`.

### Changed

- **An orifice bore against the momentum chamber it opens into (#1016).** #788's rule, applied to an adjacent `MomentumChamberNode` with areas: the bore's own geometric area (a plate's total hole area for the effusion, film, slot and impingement elements; never Cd A).
  - **Before:** a bore wider than the chamber was accepted silently. The solve then found a root where the orifice created total pressure (refused by the second-law check), left the chamber with no subsonic state, or stalled. Of the stress `random` family's 200 draws, 13 of the 19 with such a pair failed.
  - **Now:** a wider bore is refused at `resolve_topology` with a message naming the element, the chamber and both areas. A bore equal to the chamber's area (to `ORIFICE_BORE_PIPE_RTOL` on the equivalent diameter) is no orifice: zero loss, a `Pt_in - Pt_out` row. A smaller bore is unchanged. Only a chamber area that states the geometry is compared. The `DEFAULT_CHAMBER_AREA` placeholder is not, nor is an area a merge chamber inherits from an orifice main inlet, nor (with `in_pipe=False`) an area inherited from a channel.
  - **GUI:** the Inspector shows the solver's verdict for an orifice, effusion plate, film row or slot next to a momentum chamber, the pass-through or the refusal, from the resolve preview.

- **An element's explicit area against the pipe or chamber it adjoins (#1023).** The rule of #788 (bore vs pipe), for the other area-bearing elements. Equal means within `NEIGHBOUR_AREA_RTOL` (1e-6, #788's band) on the equivalent diameter.
  - **Equal is a pass-through.** An `AreaChangeElement` with F0 == F1 was refused ("use a LosslessConnectionElement"). It is now no area change: zero loss, the row `Pt_in - Pt_out`, the constant Jacobian of a lossless connection. The resolve preview reports `pass_through = 1`. A loop of such pass-throughs (or of #788's pass-through orifices) is refused at validate, as a loop of lossless connections is.
  - **Refused (BREAKING):** an explicit `AreaChangeElement` F0/F1, or a `TeeJunctionElement` F_C/F_branch, that contradicts the single pipe attached at its face. The `ValueError` names both. Before, the loss was computed for an area ratio the network does not contain. Also refused: an `EjectorElement` `mixing_area` larger than its outlet pipe, which was a contraction with no loss at the diffuser exit. A wider outlet pipe stays accepted.
  - **Warned:** a `MomentumChamberNode` smaller than an adjoining pipe is a contraction or expansion with no loss. `FlowNetwork.validate` issues a `UserWarning` naming the chamber, the pipe and both areas, and suggests an explicit `AreaChangeElement`. Never for the `DEFAULT_CHAMBER_AREA` placeholder.
  - **Unchanged:** a tee branch wider than the main (psi < 1).
  - **GUI:** the Inspector shows the resolve preview's verdict: an area change's pass-through or refusal, an ejector's refused mixing area, and a momentum chamber's warning.
  - **Stress harness:** the `ejector` family fits a mixing duct drawn wider than its outlet pipe to the pipe, with no new random draw. Of seeds 0:100, ejector 1, 11, 68 and 93 change.

- **A restriction with zero area is closed (#1024).** An `OrificeElement` with bore 0 (`diameter=0` or `area=0`), an `EffectiveAreaConnectionElement` with `effective_area=0` and a `DiameterDischargeCoefficientConnectionElement` with `diameter=0` are walls: the row is `m_dot = 0` with the constant Jacobian 1, and no kernel is evaluated.
  - **Before:** they passed `validate()` and the solve stopped at the starting point with "nozzle_flow: A_eff must be positive".
  - **Now:** they solve with no flow through them. The diagnostics report `closed = 1`, and the GUI's resolve preview reports a `closed` field. The start guess seeds a closed element's flow at zero.
  - **Topology:** a part of the network joined to the rest only by closed elements keeps solving when it has a `PressureBoundary` of its own. A part they cut off from every `PressureBoundary` is refused at `FlowNetwork.validate()`, with the closed elements named, because nothing sets the pressure of the gas trapped there. The #1004 fraction-loss check and an orifice's dead-end test treat a closed element as a wall.
  - A negative or non-finite bore, area or effective area is refused at construction. In the GUI an orifice diameter of 0 is no longer flagged before Solve.

- **Geometry and input checks refuse what cannot be built (#1024).** Each refusal names the element or node, the input and the limit. Before this change, each of these inputs was accepted. Some then failed mid-solve with a kernel message, and others solved as something else.
  - **Holes leave material between them** (`validate()`). `EffusionPlateElement` needs `pitch_y > d` and `pitch_x > d / sin(angle_deg)`, the hole's footprint on the surface. `FilmCoolingElement` needs `lateral_pitch > d`, and with more than one row `row_pitch > d / sin(angle_deg)`. `ImpingementPlateElement` needs `pi/4 / (xn_d yn_d) < 1` and `yn_d > 1`. Between rows it needs `xn_d > 1` when inline. When `pattern="staggered"` it uses the true nearest neighbour, half a spanwise pitch over: `min(sqrt(xn_d^2 + (yn_d/2)^2), 2 xn_d) > 1`. Equal is refused too. Effusion and film rows have no stated arrangement, so they are taken inline.
  - **Ducts and nodes need an area** (`validate()`). A `ChannelElement` diameter and a given Dh, and a given `MomentumChamberNode` or `CombustorNode` area or Dh, must be positive and finite. A duct of no area was refused only mid-solve ("fanno: L and D must be positive"), and a chamber of area 0 or -1 reported Mach and Re as 0. To shut a path, use a closed (zero-bore) restriction.
  - **Hydraulic diameter: Dh <= 2 sqrt(A/pi)** (`validate()`), the isoperimetric bound, for a chamber or combustor whose Dh is given (`ISOPERIMETRIC_DH_RTOL`, 1e-3 relative, absorbs an area typed to three figures). The message gives both numbers. A `ChannelElement` is not bounded: a channel that stands for a share of a wider duct, with only its part of the wall wetted, keeps the duct's Dh above the bound (the C3X example's hole segments). Its Dh must be positive.
  - **`CompressorElement` / `TurbineElement`: `pressure_ratio` must be finite.** `inf` passed the `>= 1` test. The GUI schema refuses it too.
  - **`CombustorNode(method=...)` is checked against `"complete"` / `"equilibrium"`**, at construction and at `validate()`. A typo such as `"equilibirum"` ran as `"complete"`.
  - **`EnergyBoundary(fraction=...)` must be above -1**, per boundary and summed over a node's boundaries. At -1 the node landed on the 298.15 K reference whatever it was fed, and on a combustor the solve raised a bare `ZeroDivisionError`. Below -1 the node went below the reference temperature, and from a 1000 K inflow below 0 K, where the solve failed with "thermo_state: temperature must be positive". To set an outlet temperature, use `Tt_target=...`.
  - **GUI:** the resolve preview carries each item's `refusal`, the message its `validate()` would raise, and the Inspector shows it before Solve.

- `sympy` is no longer a runtime dependency of `combaero`: `pip install combaero` needs only numpy and scipy. Nothing in the package used it at runtime since the MPCE Jacobian moved to C++ (#271). Its one importer, the sympy cross-check of that Jacobian, moved from `combaero.network._mpce_jacobian` to the test suite (`python/tests/_mpce_jacobian.py`). `sympy` is now in the `dev` dependency group, with the derivation scripts under `scripts/` that use it. A test now imports every module of the package in a fresh interpreter and checks that neither sympy nor mpmath is loaded.

### Fixed

- **Duct Mach below the 1 Pa drive floor (#1017).** `fanno_channel_flow` and `fanno_rayleigh_channel_flow` blend the flux below `FANNO_FLOW_DRIVE_FLOOR`, but returned the floor's `M_in` and `M_exit`.
  - **Before:** at a 0.035 Pa drive on 20.8 bar the reported Mach was 23x the flow's (80x at 0.01 Pa on 1 bar). A `ChannelElement`'s `mach_in` / `mach_out`, `M_in_duct` / `M_exit_duct`, its end states (`rho_ref` 1.05 against 24.2 kg/m^3) and the second-law heat band all read it. The heat band then refused a converged root (stress random 99).
  - **Now:** `M_in` comes from the isentropic inlet at the blended flux, and `M_exit` keeps the floor's exit-to-inlet ratio. `G` and its derivatives are unchanged, so residuals and Jacobians are bit-identical. The new field `FannoChannelFlow.M_exit_root` carries the solved root's exit Mach, which is the floor's when below it. `ChannelElement` warm-starts from it.

- Retry chain at a fixed budget (#1018): with the boundary-condition ramp to follow (the GUI, `NetworkRunner`, the stress harness), a capped Newton homotopy now gets at most `NEWTON_HOMOTOPY_BUDGET_BEFORE_RAMP` = 3 s instead of its 5 s, so the ramp, which is what solves those draws, starts in time; the homotopy's own stress rescues take under 2 s. A network whose preferred outlet-referenced warm start (#706) failed no longer solves from that same seed a second time in the auto-retry; it tries the inlet-referenced seed instead. Stress random 0:200 at 20 s, serial: random 119 now converges, no draw is lost.

- `fanno_channel_flow`'s choked derivatives are now exact for very short ducts. This fixes the compressible `CrossflowSegmentElement` station Jacobian at near-zero segment length (#1019). A choked duct a few microns long carries a flux within ~1e-6 of the sonic one. The choked branch differenced L*(G; Pt0, Tt0, Y) through a re-solved inlet Mach, and its steps crossed sonic. dG/dTt0 came out 17% high at L = 1 um, so the station rows' columns that move the station temperature (the segments' flows and a hot injection's flow) were 7-17% off. The branch now takes the inlet-Mach limit exactly (Leibniz), as the unchoked branch already did. Flows and residuals are unchanged; only the Jacobian moves.

- **The Levenberg-Marquardt fallback no longer turns a nearly singular Jacobian into NaN iterates (#1020).** MINPACK's LM detects a singular Jacobian only when a pivot is exactly zero.
  - **Before:** at an iterate far from the root, an effusion plate's McGreehan-Schotsch crossflow Cd underflowed to ~1e-308. The plate's row was the only one that set a chamber's pressure level, so the pivot was 1e-308 instead of 0. The step overflowed, the trust radius became NaN, and every later LM iterate was NaN. Stress effusion 34 failed with no_progress at any budget.
  - **Now:** the Jacobian given to LM has its entries below machine epsilon times their column's largest set to zero. QR cannot resolve entries that small anyway. LM then takes its rank-deficient step, and effusion 34 converges cold. hybr's Jacobian and every residual are unchanged.

- A cold `NetworkSolver.solve` of a network with an `ImpingementPlateElement` under `cd_crossflow='florschuetz_metzger_su_1982'` (CR-3630 Table 4.1) now runs the constant-Cd stage first: every crossflow-Cd plate at its own Cd, then the exact problem from that root, and the plain cold chain only if that fails (#1021). The cold seed put the stagnant rows' jets backwards in a basin too narrow for the primary attempt: CR-3630 B(5,4,1) at mc/mj 0.99 went through the failed primary and outlet-referenced attempts before the constant-Cd retry found the root (25.5 s); it now solves in 5.5 s to the same root (7e-8 relative). Cases that solved cold before reach the same root (below 1e-7 relative) but take the two stages (B(5,4,1) mc/mj 0 and 0.51: about 1.7 s before, 4 s now). Networks under the default Rogers-Hersh law, or `'none'`, keep their path unchanged.

- **A node's wetted surface area is never its flow area (#1022).** When a `MomentumChamberNode` or `CombustorNode` inherited its flow area, `resolve_topology` wrote the cross-section into `surface.area`, which is the A in a wall's h A.
  - **Before:** a given `ConvectiveSurface(area=0.7)` on a node fed by a 0.05 m pipe became 0.00196 m^2, giving 356x less wall heat with no warning. An unset surface (0) became the cross-section, so every node with an inherited area gained a convective surface. The same write was in the tee and MPCE port sizing, the ejector's outlet sizing and `PressureLossElement`'s area inheritance.
  - **Now:** resolve never writes `surface.area`. A given wetted area is kept and an unset one stays 0 (no convective surface), so a `ThermalWall` on such a node is refused by name. Set the wetted area to put a wall on it.
  - **Effusion gas side:** an `EffusionPlateElement` discharging into a chamber reads the chamber's h whether or not the chamber has a wetted area of its own. The plate owns the wall, and its wall temperatures used to rely on the inherited cross-section.
  - **GUI:** combustors, momentum chambers and discrete losses have no wetted-area field. Their surface used to be the flow area. Now a wall drawn to one gives that side its area (the side's `area_a`/`area_b`, else the wall's area). The coupling already used that area, so wall results do not change. An unwalled one no longer reports h, Nu or T_aw.

- A momentum chamber or combustor whose heat transfer is not evaluated (no convective surface, no flow, or a choked-bounded junction) reports `Nu`, `St`, `htc` and `T_aw` as NaN, which the GUI shows as no value. Since 0.11's #1022 this covers every chamber and combustor without a wall; before, they read 0 (zero heat transfer) and `T_aw` = Tt.
- GUI: an orifice that is a pass-through against its chamber (#1016) no longer shows beta = 1 and E = 1 next to the "beta = 0" hint; beta and E show only against a pipe diameter.

## [0.10.0] - 2026-11-06

### Added

- **A fraction loss without `m_ref` on its zero-flow blend is flagged (#1004).** `ConstantFractionLoss` and `LinearThetaFractionLoss` without `m_ref` hold dPt = xi Pt at every flow but zero. Where the network gives one a smaller drop, its only root is the forward/reverse blend at |m_dot| < `FLOW_REVERSAL_M_SCALE`, and the solve converged there reporting success with nothing to say the loss law is not met.
  - **Now:** every `PressureLossElement` reports `fraction_unanchored`: 1 for such a loss on its blend when its flow component is driven and the loss is not a dead end, 0 otherwise, NaN where its flow is not finite. `__fraction_unanchored__` lists them. Reported only, never a demotion.
  - **Driven, from the inputs only:** pressure boundaries whose Pt differ by more than 1e-9 of the largest, a `MassFlowBoundary` feeding a nonzero flow, a `FuelBoundary` whose target is above the coldest other boundary Tt, or a compressor or turbine at a pressure ratio other than 1. **Dead end, from the topology:** one side of the loss reaches no pressure boundary, fuel controller or nonzero feed without it. Checked once per solve, after every retry has restored the boundaries.
  - **A failed solve** names the losses flagged at the best iterate in `__message__`.
  - **GUI:** the discrete loss's Inspector has a "Use the solved flow as ṁ_ref" button. It is enabled after a converged solve in which the loss carries flow. At ṁ = ṁ_ref the loss is exactly ξ, so the solved point is unchanged. A flagged loss shows a note pointing at the button. Nothing changes by default: ṁ_ref stays unset and the palette is unchanged.

### Changed

- `FlowNetwork.validate()` now refuses a path between two `PressureBoundary`s made only of elements whose total-pressure relation does not depend on their own flow -- the edge set of the closed-loop check (lossless connections, free vortices, a splitter's reference outlet, fixed-ratio compressors/turbines, fraction losses without `m_ref`) (#1003). Each boundary fixes its Pt, so such a path's pressure ratio is fixed at both ends: in the GUI, a plenum with plain edges to two outlets at equal Pt "converged" with success at arbitrary splits, some with an outlet flowing backwards, and at Pt 0.1 % apart it had no root and failed at about 6e11 kg/s through the pair. The message names the two boundaries and each element on the path, and says how to fix it: put a flow-dependent loss (a duct, an orifice, a head loss) on each outlet path, or merge the outlets into one boundary. A single such link into one boundary is unaffected.

- **An effusion plate with no coolant flow is flagged, never silent (#998).** Since #995 its wall keys were NaN on a plenum or a stagnant chamber and the solve still reported success; in hot gas the adiabatic wall T came with no warning that it melts the plate.
  - **No cooling path refused.** No coolant through the holes and no gas flow on the other side (a plenum, or a chamber with none): the wall T is undetermined, so `no_cooling_path` is 1 and a converged solve is demoted to INCONSISTENT with `EffusionPlateElement 'id': no cooling path: wall T undetermined (heats until failure)` in `__message__`, the plates in `__wall_undetermined__`, as a `Tt_target` with no through-flow is (#792).
  - **Temperature limit.** `EffusionPlateElement(material=...)` / `EffusionLiner(material=...)` name the plate's alloy for its limit only: `T_wall_limit` is the hot end of its #882 envelope (`cb.GENERIC_METAL_SOLIDUS`, 1533.15 K, without one, now exposed), and `wall_above_limit` flags a face above it (NaN when the wall T is undefined), listed in `__wall_above_limit__`. The conduction stays `wall_conductivity`. The GUI plate and liner take the material and show both notes.
  - **No flow below the floor.** The plate counts as having no coolant, and its gas-side chamber as having no gas flow, below `cb.MIXER_HEAT_MDOT_FLOOR` (1 mg/s), not only at exactly 0: a zero-flow root ends at ~1e-31 kg/s through the holes, where the coefficient put the cold face at 5.7e11 K, and ~1e-13 kg/s in the chamber, which made a stagnant liner an adiabatic 1500 K wall.
  - Not modelled, by decision: radiation, backside convection without hole through-flow, conduction to the mounts. The code comment at `wall_heat_transfer` says how each could be; #898's wall-to-wall conduction is the path for a locally starved plate.

- `FlowNetwork.validate()` refuses a closed loop made only of `LosslessConnectionElement`s -- two in parallel between the same nodes, or a ring of them -- and names its elements (#988). Total pressure is equal all round such a loop, so any circulating flow is a root: the solve returned an arbitrary one, start- and release-dependent (stress random 5: 159 kg/s circulating while the network carried 8.7 kg/s), and a heated element on the loop got an arbitrary heat transfer. In the GUI this is two plain edges drawn back and forth between the same nodes, or a ring of plain edges; the solve now says so. A loop with any loss element in it is determined and unaffected. Remove one of the connections, or give one path a loss element.

- **The cold start no longer shares a temperature across wall-joined flow components (#899).** Each flow component's reference state now comes from its own boundaries alone, its reference Tt included; before, components joined by thermal walls started at their thermal group's mean Tt. Since #990 moved the resistance seed to each component's own Tt, the group mean reached only the propagated start and the reference m_dot.
  - Networks of one flow component, or of components no wall joins, start bit-identically.
  - Wall-joined networks start slightly differently (stress: hot/cold pairs by at most 1e-13, impingement median 1e-7, pin fins by up to 46 % on the propagated start) and converge to the same roots (564/564 stress seeds, all 62 NASA TN D-7156 runs, C3X), at the same wall-clock time within 1 %.

- `FlowNetwork.validate()` now refuses any closed loop made only of elements whose total-pressure relation does not depend on their own flow, not just a loop of lossless connections (#988): `LosslessConnectionElement`, `VortexElement` (a free vortex conserves Pt), a flow splitter's reference outlet, `CompressorElement`/`TurbineElement` at their fixed pressure ratio, and a `PressureLossElement` with `ConstantFractionLoss` or `LinearThetaFractionLoss` and no `m_ref`. Round such a loop the pressure ratios either fail to close, so there is no solution (a fraction loss in parallel with a lossless link "converged" only by sitting in its zero-flow blend, or failed with a circulation of thousands of kg/s), or close exactly, so the split between the paths is not determined. The message names each element and why, and says how to fix it: put a flow-dependent loss (orifice, channel, head loss) on one path, or give the fraction losses an `m_ref`. A fraction loss with `m_ref` fades to zero at zero flow and is not affected. A loop made only of lossless connections keeps its #988 message.

- A result's keys no longer depend on the operating point: a value that is undefined at a point is reported as NaN (Python) / `null` (GUI `/solve`) instead of the key being left out. This covers the second-law band (`sigma_norm`, `sigma_norm_low`, `sigma_norm_high`, now NaN on a `not_checked` verdict, e.g. at zero flow), a channel's `Re_surface` / `f_surface` / `surface_extrapolated` and a discrete loss's `Nu` / `St` / `htc` / `T_aw` / `f` without flow, the ejector's `P_c_star` / `P_b0` / `entrainment_ratio_critical` unchoked and `head_ratio_N` / `efficiency` without a drive, the tees' `port_i_K`, `port_i_flow_fraction`, `K_straight` / `K_branch` / `K11` / `K12` and `mass_flow_ratio` in every flow pattern (ideal momentum, Mynard, constant K), and the effusion plate's coolant- and gas-side keys without coolant flow. Code that tested `key in result` to detect a regime should test `math.isnan(result[key])`.
- `jet_flow_nonuniformity` (impingement array) and `bleed_nonuniformity` (effusion liner) are NaN when no row flows or one runs backwards (they were inf); inf remains for a starved row beside a flowing one.
- `scripts/gui_result_keys.py` counts a key whose value is `null` as present: the generated list is the result schema, not the values one solve produced.

### Fixed

- **An orifice in a pipe solves at zero flow (#997).** With an approach pipe (beta > 0), the Jobson model (the default for a sharp-edged bore, and the Lichtarowicz, IdelchikThick and McGreehanSchotsch correlations) failed at equal boundary pressures with `chamber_static_state: no convergence ... (last step [K] = nan)`.
  - **Cause:** two roundings. `duct_static_state` (the approach pipe's static state) came back one ulp above its own total pressure at zero flux. The Jobson momentum balance wrote the jet's `1 - (p/p0)^((n-1)/n)` as a difference near 1, which rounded to zero a few ulp below `p0`. At zero flow the kernel then divided a positive drive by a zero jet velocity, and the infinite flow it returned reached the next static-state solve as a NaN.
  - **Now:** the static pressure is solved relative to `Pt`, so it equals `Pt` to the bit at zero flux and never exceeds it. The jet term is evaluated as `-expm1(((n-1)/n) log1p((p - p0)/p0))`, so the flow follows `m ~ sqrt(dp)` down to a one-ulp drive. Solves away from zero flow change only at the rounding level.

- **The cold start's resistance seed runs at each stream's own temperature (#990).** In a network of flow components joined by thermal walls, the resistance-network seed now uses the component's own mean boundary Tt, for both its density and its compressible kernels. Before, it used the thermal group's mean (#899). Every other start estimate still uses the group's mean.
  - **Before:** NASA TN D-7156 run 19 seeded its Mach 0.7, 297 K coolant passage at the group's 446 K. Its pressures started 19 % high (172.8 against a 145.2 kPa root) and the default cold solve needed the warm-start retry.
  - **Now:** all 62 D-7156 runs converge cold, and run 19 takes about a third of the time. Only networks of several flow components joined by walls start differently, and their solved roots are unchanged.

- **The cold start seeds a fraction pressure loss from its own law (#989).** In the resistance-network seed, a `PressureLossElement` with a `ConstantFractionLoss` (or a `LinearThetaFractionLoss` with no theta source) now takes at most its own drop, xi Pt_up. Once its reference-flow resistance would take more, it enters as the pressure ratio Pt_dn = Pt_up (1 - xi) and leaves its flow to the rest of the network. Below a design flow `m_ref` it enters as the dynamic head the loss fades with.
  - **Before:** on stress random draw 3832186739 the seed sent 0.56 kg/s through a 1 % fraction loss that carries 117 kg/s at the root, and the cold solve failed.
  - **Now:** that draw starts within 1.5x of its root flows and converges cold. Only networks with a fraction loss start differently, and their solved roots are unchanged.

- GUI: a fraction loss's ṁ_ref set from the solved flow is shown at 6 significant figures; the stored value stays exact, so the solved point does not move (#1006).
- GUI: element cards show a mass flow below 1e-3 kg/s in exponential notation, so a blocked fraction loss at 5.5e-5 kg/s no longer reads "0.0001 kg/s".
- A converged solve whose root is refused (outcome inconsistent, e.g. an effusion plate with no cooling path) no longer gets "[... auto-retry skipped: neither ... proxy solve converged]" appended to its message, and the GUI's Live Telemetry and probes say the root was refused instead of "did not converge". `/solve` now returns `converged`, `consistent` and `outcome`.
- GUI: `T_wall_limit`, `wall_above_limit`, `no_cooling_path` and `fraction_unanchored` have readable labels in Live Telemetry and Display Fields.

- A node with an energy boundary now always reports `{node}.Q_withheld`, 0.0 when nothing is withheld. Before, the key appeared only when the value was nonzero, and `m * (1/m)` rounding to `1 - 1 ulp` (about one flow in eight, e.g. 0.09 or 0.36 kg/s) reported ~1e-11 W of "withheld" heat on a flowing node. Which keys a result has, and the GUI's result-key list, therefore depended on the last bit of a solved flow and differed between platforms. The rounding remainder (below 1e-12 of the heat given) is now reported as 0.0; heat withheld by a stagnant node is unchanged.
- `{node}.Q_fraction` is likewise always reported on a node with an energy boundary: 0.0 when there is no fraction boundary or no flow for it to be a fraction of (it was absent).

- An `EffusionPlateElement` on a gas-side chamber with no coolant flow reports the adiabatic limit (q = 0, both faces at the gas side's adiabatic wall temperature, no film) instead of dropping its wall keys; on a plenum, or on a chamber with no gas flow, the wall temperatures are NaN. A chamber with no gas flow no longer switches the plate to the plenum key set.
- The GUI's flow splitter always reports `valve_dPt_max` and `fraction_error_max`, NaN (null) when no outlet row carries one (no flow to split); `fraction_error_max` was absent at zero flow on Linux.

## [0.9.0] - 2026-11-04

### Highlights

The 0.9 release (66 commits since 0.8.0), the first from the `thirox-lab`
organisation; the sections below are the record. What changes numbers or code
on upgrade:

- **The 0.8 deprecations are removed.** The old stagnation spellings
  (`T0`/`P0`/`T_total`/... and the `*_from_static` names) now raise; use
  `Tt`/`Pt`. The `"<chamber>.T_static"` key, Stolz (#876) and
  `k_aluminum_6061` are gone.
- **Invalid input raises `ValueError` everywhere (#928)**, not `RuntimeError`,
  and NaN input is refused rather than answered (#769).
- **Materials (#881, #882).** Outside a material's sourced range k(T) follows
  an envelope instead of holding the end value, flagged by
  `material_in_range`; `ss316` is the Gen3 CSP measured table; `k_haynes230`
  is the brochure's alpha * cp * rho.
- **Inner solves reach their roots (#868, #927, #935-#943)**: cap exits are
  reported instead of returned, so results move slightly.

New: compressor and turbine elements, a fuel controller, a flow splitter,
liner radiation, ThermalWall per-side areas and a side with no stream
(adiabatic or heat loss), `ss304`/`ss310` materials, the NASA TN D-7156 and
C3X vane cases, GUI nodes for the new elements, and cold starts that converge
without the retry on channel chains (#971).

### Added

- `ThermalWall` takes a side with no stream (#970): `element_b=None` with
  `side_b="adiabatic"` (no heat; the wall sits at side A's T_aw, or with
  `GasRadiation` on side A at the face's radiation balance, C++
  `wall_adiabatic_face`) or `side_b=HeatLoss(T_ext=..., h=...)` (side A's
  stream exchanges with a fixed external temperature through the series
  resistance, epsilon-NTU against a reservoir, analytic Jacobian). The heat
  lost is reported as `{wall}.Q_loss` and booked by the ledger
  (`__ledger__["Q_loss"]`) and the wall's second-law check. The NASA TN
  D-7156 runner (#764) models its stagnant backside with it instead of a
  token flow. Not in the GUI yet.

- **`calc_T_from_cp_roots`: every root of cp(T) = target (#963).** cp(T) is not monotone (N2 has a tail maximum at 127 K and a minimum at 215 K; O2, H2O, CH4, CO and NH3 have tail minima near 100-125 K, CO2 at 6.2 K). `calc_T_from_cp` still returns the root nearest `T_guess`, unchanged bit for bit. The new opt-in call needs no guess.
  - **Returns:** a sorted list of `CpRoot` (`T` [K], `dT_dcp` = 1/(dcp/dT), `double_root`) in [`T_min`, `T_max`] (default [1, 6000] K). The list is empty when there is no root.
  - **Method:** cp is split at its extrema: the sign changes of the analytic dcp/dT on the scan nodes, bisected, plus a pair inside one step, which a search finds where |dcp/dT| dips. Each monotone piece then holds at most one root. This finds the close pairs `calc_T_from_cp` reports as "no root" (O2's cp(127 K) is reached again at 120.8 K).
  - **Double root:** a target equal to cp at an extremum (to 1e-13) is returned once, with `double_root=True` and `dT_dcp` NaN.
  - **Errors:** `ValueError` for invalid input or range, and for argon's exactly constant cp, whose roots are not isolated.

- **Film-side heat transfer coefficient near a slot, opt-in (#962).** `SlotFilmElement(..., film_htc="seban1960")` gives the wall side its film covers the film's own h, from the slot flow and the distance from the lip, in place of the duct's. T_aw still comes from the film effectiveness.
  - **Seban (1960) Fig. 12 as drawn**, on slot-fluid properties: `St_s = 0.135 Re_s^-0.3 (x/s)^-0.3` near the slot and `0.40 Re_s^-0.3 (x/s)^-0.6` beyond, C1-blended at their x/s-37 crossing. Seban's Eq. 13 prints 0.10; his own curve B, Fig. 9 and data (digitised, 0.134) contradict it. Lefebvre & Ballal's Eq. 8.41 (0.069) carries that misprint forward.
  - **Element mean** over the wall's span, not a point value. The side hands back to the duct's h past x/s 300 and below m = 1 (Seban's flat-plate finding; C1 over 0.1 in m). Outside his data otherwise (x/s < 4, Re_s 1.4e3-8e3, m up to 8.91) it is applied and flagged.
  - **Named option** `"ballal_lefebvre1973"`: the book's m > 1.3 wall jet `Nu_x = 0.10 Re_s^0.8 (x/s)^0.44` (Eq. 8.44). It reproduces the book's worked example but sits 20-67 % above the Cranfield rig's own data, so it is not the recommended form.
  - **Kernels:** `slot_film_stanton_seban_1960`, `slot_film_nusselt_ballal_lefebvre_1973`, `slot_film_htc_seban_1960`, `slot_film_htc_ballal_lefebvre_1973` (analytic partials in G_s, T_s, M, h_duct); constants `SEBAN1960_C_NEAR` / `_FAR` and their `_PRINTED` values, `SLOT_FILM_BRANCH_*`. Analytic Jacobian through the slot flow, coolant temperature and gas flow. Reported: `{wall}.film_h`, `film_h_duct`, `film_h_film`, `film_h_branch`, `film_h_in_box` and parts. In the GUI's slot-film inspector too.
  - **Scored:**
    - Fidelity on Seban's own Fig. 12 points (38, 5 clean): bias -1.0 %, RMS 11 %.
    - Cross-source accuracy on Cranfield SME4 Fig. 1 (Hartnett; Pai & Whitelaw): bias -2.5 %, RMS 18 %.
    - The book's film-cooled and TBC liner examples as networks with #960's radiation: 1285/1267 and 1308/1240/1224 K with Seban, 1269/1252 and 1288/1224/1209 K with the wall jet, against the printed 1283/1265 and 1304/1236/1220 K. This is a method comparison, not fidelity.
  - **Default off:** without `film_htc` every network is bit-identical. The digitised Seban Fig. 12 and SME4 Fig. 1 points are in `validation/cooling/data/`.

- **`FlowSplitter`, a flow splitter at given fractions (#903).** One inlet node to N outlet nodes, each outlet carrying a given fraction f_i of the throughput: a mass-flow controller per outlet. Changing a network's bulk supply is one number, not one `MassFlowBoundary` per path.
  - **Fractions:** each in (0, 1), summing to 1 within `SPLITTER_FRACTION_SUM_TOL` (1e-9). A set that does not is refused, never renormalised. One outlet may be the remainder (`None`, 1 minus the others).
  - **How it solves:** one `SplitterOutletElement` per outlet, a two-port from the inlet node with one unknown, its `m_dot`. A metered outlet's row is m_i - f_i sum_j m_j = 0, linear with an exact Jacobian. It holds no pressure, so the outlet's total-pressure drop is free: the drop its valve would have to take.
  - **The reference outlet:** one outlet has no valve (Pt_out = Pt_in), since only N - 1 fractions are independent. It must be the path that needs the most pressure, and the solve finds it (`auto_reference`, default): it starts from the given `reference`, else the remainder, else a first guess. When a metered outlet would need a pump, the solve moves the reference to the outlet needing the most pressure and solves again. This is an outer loop of at most N - 1 swaps, reported in `__reference_swaps__` and the message. With no feasible reference, or `auto_reference=False`, the solve is refused, and the message names the outlet to make the reference.
  - **Reported per outlet:** `fraction`, `fraction_achieved`, `valve_dPt`, the metering effective area `valve_CdA` (an isentropic jet from the inlet's stagnation state, `valve_choked`), and with a reference area `valve_Cd` and `valve_zeta` (the K). `FlowSplitter.summarize` collects them.
  - **Thermo:** every outlet carries the inlet's Tt and composition; Pt follows each outlet's own path.
  - **Tested:** fractions to 1e-13; drops and metering against a hand calculation; every Jacobian row against central differences. The C3X example's 19 coolant paths through one splitter reproduce its per-path flows and the pressures its per-path boundaries float to. A many-inlet `PlenumNode` and merge `MomentumChamberNode` (300-1900 K, four compositions) close a hand energy and species balance.

- **`gui.backend.export.network_to_gui(net)`: a `FlowNetwork` built in a script becomes a GUI network file (#902).** It opens in the GUI and solves to the same results.
  - Supported: boundaries, plenum, momentum chamber (side streams on its side inlet), combustor, channel, orifice, pressure loss (all four correlations), lossless connection, thermal wall, smooth surfaces.
  - Layout: by flow depth, given positions, or a callable.
  - Compositions are written as exact mass fractions.
  - Anything without a GUI form is refused with `UnsupportedForGUIError`, which lists it.
  - The C3X worked example now writes its GUI file through it, in place of about 150 lines of hand-built nodes and edges.

- **Network results report the Stanton number `St`** for every convective surface: channels, impingement plates, pressure-loss elements and chamber or combustor nodes with a surface (#895). It is St = Nu / (Re Pr) at the heat-transfer correlation's own Re and Pr: the jet Re for impingement, the pin Re for pin fins. That makes it the h the element reports in dimensionless form, which is how turbine heat transfer is usually quoted. It is 0 where no Re Pr is defined. Diagnostics only: no residual changes.

- **`cb.compression` and `cb.expansion` (#893).** They give the outlet total temperature, the work per unit mass (`dh`) and both efficiencies of a compressor or turbine path at a total-to-total pressure ratio. Pass exactly one of `eta_poly` or `eta_is`. The analytic Jacobian of T2 and dh in T1, PR, the efficiency and the mass fractions is computed in C++.
  - **The polytropic form is exact for a variable cp:** S(T2) − S(T1) = (R/η_p) ln PR, with no path integration.
  - **Checked against Cantera:** isentropic to 3e-8 (humid air at PR 2–40, combustion products at PR 2–20); polytropic to 1.6e-9 against the path ODE on Cantera's cp. Argon, whose cp is 5/2 R, reproduces the constant-gamma closed forms.
  - **Refused:** PR < 1, an efficiency outside (0, 1], and neither or both efficiencies given.
  - These are the kernels for the compressor and turbine network elements that follow.

- **`FuelBoundary`, a fuel controller (#892).** It finds the fuel flow that puts a combustor's exit at a target total temperature, replacing an outer loop around the solve. The fuel flow is a solver unknown (as ln m_dot, so it stays positive), closed by one temperature row, Tt(combustor) - Tt_target, with an analytic Jacobian through the combustor's own Tt relay.
  - **Route:** `combustor` names the controlled combustor, like a probe. The fuel may reach it through any elements -- an injector, a lossy fuel line, a premix plenum the air also enters -- and gives the direct-injection answer to 1e-6. Refused before the solve: a route that lets part of the fuel bypass the combustor, a second controller or a `Tt_target` heat boundary on the same combustor.
  - **Branch:** a reachable target has a lean and a rich root. `branch="lean"` is the default; `"rich"` is for RQL staging or reformers. The start is put on that branch (complete combustion). A solution on the other side of stoichiometric is refused with the reason, never returned.
  - **Refused, with the reason:** a fuel composition with nothing combustible (N2, CO2, H2O, Ar) or an oxidizer (O2), before the solve; a target not above what the combustor reaches with no injected fuel (the oxidizer's Tt, or its own burn when it carries fuel); a target above the peak, reported at stoichiometric.
  - **Reported:** `<fuel>.m_dot`, `<fuel>.phi` (the reactants' combustibles over their O2, so humid air's water vapour does not count as fuel), `<fuel>.Tt_target`, `<fuel>.Tt_without_fuel`, `<fuel>.rich_branch`.
  - **Checked on every pure species:** the ten fuels (the seven C1-C7 alkanes, CO, H2, NH3) meet lean and rich targets that the combustor's own model reproduces to 1e-3 K. Air that carries CH4 of its own works too.

- **`ThermalWall` per-side areas (#887).** `area_a` / `area_b` give a wall whose two faces differ, such as a hole in a thick wall, a finned face or a curved wall. Each side's convective resistance is formed on its own face: 1/UA = 1/(h_a A_a) + Σ t/(k A_m) + R_fouling/A_a + 1/(h_b A_b).
  - **`wall_area`** chooses the layers' conduction area. `"log_mean"` (the default) runs the area linearly through the stack, with each layer on the log mean of its two faces; this is exact for a thick cylinder, layer by layer. `"slab_a"` and `"slab_b"` put a plain slab on one face.
  - **Refused:** both areas or neither, and not together with `contact_area`. A wall without them is unchanged.
  - **Checked:** against the analytic thick cylinder to 1e-13, single and two-layer. The Jacobian matches finite differences, and equal faces reproduce the one-area wall.
  - **GUI:** the thermal edge's "Per-side areas" option (side A, side B, conduction area). `network_to_gui` exports them.
- **The C3X example's walls now carry their two faces** (the segment's on the gas side, the hole's share of its perimeter on the coolant side). The coolant channels' `Nu_multiplier` area factor is gone. Tw is the same to 3e-11 K.

- **`ss310` material and `k_stainless_steel_310` (#881)**, ASTM 310 / 1.4845 heat-resisting stainless steel. Source: Georgia Tech Gen3 CSP measured k, valid 533-1073 K (260-800 C), 95% band 0.8-2.5 W/m K. It is not joined to a room-temperature value: NIST's cryogenic fit (11.6 W/m K at 300 K) and this data's trend (about 13.7 there) are different samples.

- **A network validation case: NASA TN D-7156 (#764).** Marek and Juhasz (1973) cooled a plate in a combustor's exhaust with a slot film, a backside convection passage, or both. All 62 runs are now built as combaero networks: a `SlotFilmElement` (Seban film-side h), ten filmed `ThermalWall` segments, a backside channel chain, and vitiated hot gas.
  - **The data:** Table II in SI and US units, transcribed twice with 0 disagreements and signed off 2026-11-03 (`validation/cooling/data/marek_juhasz1973/`).
  - **The scores** are cross-source accuracy, with nothing tuned. They are staged against the report's +-11 K band. Wall bias:
    - convection only: -49 K;
    - film only: +25 K;
    - combined: +2 K, with RMS 44 K (+37 K at m >= 1, -22 K at m < 1).
  - **Also reported:** the film-side h comparison (seban1960 / none / ballal_lefebvre1973) and thermocouples 1-3 apart from 4-10.
  - **Where:** `validation/cooling/marek_juhasz_network_runner.py`; pinned in `test_marek_juhasz_d7156_network.py`.

- **Liner radiation, opt-in (#793).** Lefebvre & Ballal (2010) Ch. 8, gray and lumped, built only from sources we hold.
  - **Kernels**, each with analytic partials: `gas_emissivity_luminous` (Eq. 8.11, Reeves' Eq. 8.9 with L = 1), `luminosity_factor_hydrogen` (Eq. 8.15, L = 336/H^2), `liner_gas_radiation` (R1, Eq. 8.7), `liner_casing_radiation` and `liner_casing_exchange_factor` (R2, Eq. 8.17), `mean_beam_length` (Eq. 8.10, 3.4 V/A; `BEAM_LENGTH_COEFF_HOWELL` = 3.6 for the textbook value) and `elemental_fuel_air_ratio` (q from a mixture's C and H). Constants `STEFAN_BOLTZMANN`, `REEVES_EMISSIVITY_COEFF`, `LUMINOSITY_HYDROGEN_COEFF`, `BEAM_LENGTH_COEFF_LEFEBVRE`.
  - **On a ThermalWall:** `radiation_a` / `radiation_b` take a `GasRadiation` (gas to wall, from the side's stream Tt, static P and composition unless given) or a `CasingRadiation` (wall to a casing at the side's stream T unless given). The local wall balance is solved for the two surface temperatures and folded into each side's conductance and driving temperature, so the streams' exchange runs unchanged. Analytic Jacobian, k(T) layers included. Reported: `q_rad_a/b`, `Q_rad_a/b`, `T_surface_a/b`, `eps_g_a/b`.
  - **Default off:** a wall without a spec is bit-identical.
  - **Fidelity (printed quantities):** the book's uncooled, film-cooled and TBC liner examples are reproduced from their printed inputs (1641.6/1604.4 K against 1640/1603; 1283.0/1264.5 against 1283/1265; 1304.0/1236.6/1220.4 against 1304/1236/1220).
  - **Limits:** a kerosine-product emissivity (not for natural gas or hydrogen products), one beam length, no soot model beyond L, no spectral or wall-to-wall view factors. Beam length, area ratio and emissivities are the user's inputs. Python API only; not in the GUI yet.

- **Each material's solidus, and `get_material_info(name)` (#882).** The database records the solidus its source prints, the lower end of the melting range, together with that source. Ten materials have one: Haynes 230, 188, 263 and 25, Hastelloy X, Inconel 625 and 718, Nimonic 75, 310 and copper. ss316, x22crmov12_1, nicocraly, alumina, ysz and al6061 do not: no held document prints a value, so it is NaN.
  - `get_material_info` returns a material's sourced conductivity range, source, notes, solidus and solidus source.
  - Data only: no conductivity changes. The solidus becomes the hot end of the k(T) envelope that #882 adds next.

- **GUI fuel boundary (#892).** A Fuel Boundary node in the palette: target Tt, lean/rich branch, fuel Tt and composition (CH4 by default), an optional flow guess. A dashed yellow link shows the combustor it controls; drag its link handle onto a combustor, or connect it straight to one. Its card shows the solved fuel flow, and the Inspector shows phi and the combustor's Tt without the fuel. `network_to_gui` exports a `FuelBoundary` as this node; before, it went out as the mass-flow boundary it subclasses, which fixed the fuel flow.
- **Fuel-fixed staged start (#892).** A cold solve with a fuel controller whose air comes through a pressure boundary first holds each fuel flow at its start and solves the flows, then re-sizes the fuel on the air that root carries and solves the controllers from there. The plain cold chain runs only if that fails. Any other failed cold solve with a controller gets the same stage as a retry. The fuel start is sized on a guessed air flow, which on the C3X rig's burner was 2700 times too small: the cold solve stalled. That network now solves in about 15 s under the GUI's 30 s budget (#899 removed the fuel flow from the cold-start reference).

- **`CompressorElement` and `TurbineElement` (#893).** Network elements at a fixed total-to-total pressure ratio: Pt_out = PR·Pt_in for a compressor, Pt_in/PR for a turbine. The network sets the flow. Each takes exactly one of `eta_poly` and `eta_is`.
  - **Outlet state:** the stream reaches the downstream node at `cb.compression` / `cb.expansion`'s outlet Tt. The Jacobian chains Tt_out through the inlet Tt and composition, so a turbine fed by a burner's products is exact.
  - **Reported:** `power` [W] > 0 (absorbed by a compressor, delivered by a turbine), `dh` [J/kg] signed into the flow, `Tt_out`, `Tt_out_isentropic`, both efficiencies and `pressure_ratio`.
  - **Energy accounting:** the network energy ledger and the second-law check count the shaft work as enthalpy with no entropy (`__ledger__["W_shaft"]`). An ideal machine generates no entropy; a real one generates some.
  - **Refused, with the reason:** a root with the flow reversed (a surging compressor, a turbine driven backwards); a machine between two pressure boundaries (nothing sets its flow); a pressure ratio below 1; an efficiency outside (0, 1]; neither or both efficiencies.
  - **Tested:** a Brayton cycle (compressor, burner on a `FuelBoundary`, turbine, nozzle) closes its energy balance, and every Jacobian row matches central differences.

- **GUI compressor and turbine nodes (#893).**
  - **Inputs:** pressure ratio, efficiency type (polytropic or isentropic) and efficiency. The card shows the solved shaft power; the Inspector shows Tt in and out and both efficiencies.
  - **Defaults:** compressor PR 4, turbine PR 3.8, so a pair dropped together leaves 5 % of Pt for the combustor and ducts between them. At PR 3 (25 %) a default chain chokes its hot channel and the cold start stalls.
  - **Export:** `network_to_gui` writes `CompressorElement` and `TurbineElement` as these nodes.
  - **C3X example:** its compressor, until now computed in the script, is a `CompressorElement`, with ambient humid air at an assumed isentropic efficiency of 0.85. The scorecard is unchanged.

- **GUI: a Flow Splitter node (#903).** One inlet to N outlets at given fractions, the GUI form of `FlowSplitter`. The card grows with its outlets and shows each one's fraction, and after a solve the drop its valve takes.
  - **Inspector:** edit each fraction, pick the remainder and the starting reference (the outlet whose valve is open; the solve moves it when needed, and the card and the results show where it ended), give optional valve areas, and add or remove the last outlet.
  - **Results:** a table per outlet: flow, fraction reached, valve drop, effective area and zeta.
  - **Checks before a solve:** the fractions must sum to 1.
  - **Export:** `network_to_gui` writes a script's splitter as this node.

- **Worked example: the Hylton 1983 C3X cooled vane** (NASA CR-168015), as a script (`python/examples/papers/hylton1983_c3x_vane.py`) and a GUI network (Examples -> "Hylton 1983 C3X cooled vane").
  - The vane surface is split into 19 segments, one per nearest cooling hole. Each has a gas channel from the inlet total state to the measured local static pressure, an ss310 wall, and the hole's share of its coolant.
  - The gas-side Nu multipliers are fitted to run 158's measured wall temperatures.
  - Untuned, the smooth-channel network misses run 158 by MAE 28 K (worst -85 K at the leading edge), labelled cross-source. Run 149 with run 158's multipliers gives MAE 12 K: the multipliers absorb the run's 8.3 % turbulence, and run 149 has 6.5 %.
  - The report's h is inferred by its own FE model, so it is shown for comparison and not used as a target.
- **A third example label, "tuned"** (GUI example panel): the network carries tuners fitted to the source's own data, and its untuned error is reported beside it.

- **`ss304` material and `k_stainless_steel_304`**, AISI 304 / 304L austenitic stainless steel. Source: Ho & Chu 1977, CINDAS Report 45 (DTIC ADA129160) Table 2, recommended values with a stated uncertainty of about 5%, interpolated linearly; valid 200-1672 K, from the #882 envelope's low end to the solidus (the table's melting range, 1672-1727 K). Beyond those ends the database lookup holds the end value, and `extrapolate=True` continues the end slope. Cross-checked against NIST's cryogenic fit (-2.0% to +2.8% at 200-300 K) and the Cleveland-Cliffs 304/304L bulletin (+0.2% at 100 C, -3.7% at 500 C).

### Changed

- **Invalid input raises `ValueError` everywhere (#928). Breaking for code that catches `RuntimeError` on bad input.** 76 input checks in the C++ kernels threw `std::runtime_error` (Python `RuntimeError`) while every other module threw `std::invalid_argument` (Python `ValueError`) for the same kind of mistake.
  - **Before:** a non-positive, out-of-range or non-finite argument raised `RuntimeError` in the acoustics functions (`liner_sdof_impedance_norm`, `liner_2dof_serial_impedance_norm`, `orifice_impedance_with_flow`, `quarter_wave_resonator_tmm`, `is_whistling_risk`, `can_annular_eigenmodes`, `annular_duct_eigenmodes`, `annular_duct_modes_analytical`), the materials functions (every `k_*`, `get_material_conductivity`, `material_in_range`, `get_material_info` for an unknown name), the humid-air functions (a relative humidity outside [0, 1], a non-positive pressure, a vapour pressure above the total: `vapor_pressure`, `humidity_ratio`, `water_vapor_mole_fraction` and everything built on them, `HumidAir.set_TP_RH`, `species.humid_air_mass`), the combustion helpers (`fuel_lhv_*`, `oxygen_required_per_*`, `dryair_required_per_kg_fuel`, `equivalence_ratio`, `complete_combustion`, `set_equivalence_ratio_mole`/`_mass`, the Bilger mixture-fraction functions: wrong vector size, fractions not summing to 1, an oxidizer without O2, a fuel needing no O2), `mole_to_mass`/`mass_to_mole` (an all-zero composition), `prandtl`, `reynolds`, `peclet`, `adiabatic_wall_temperature`, `cooled_wall_heat_flux`, `nozzle_thrust` (an empty solution) and `adiabatic_T_wgs` (an empty `X_in`).
  - **Now:** they raise `ValueError` (C++ `std::invalid_argument`) with the same message. Catch `ValueError` instead of `RuntimeError` for bad input. In the GUI a refused input (a boundary's relative humidity of 1.5, for example) is now a 400 naming it instead of a 500.
  - **Unchanged:** `RuntimeError` still means a solve that failed on valid input: an inner solve's iteration cap (#868), a Fanno march out of steps, a bracket without a sign change, the dewpoint and wet-bulb iterations, `solve_orifice_mdot`'s Cd iteration, and a target no state reaches (`solve_A_eff_from_mdot`, `solve_P_back_from_mdot` above the choked flow, `solve_Pt_from_mdot`, `set_fuel_stream_for_Tad` and `set_oxidizer_stream_for_Tad` outside the achievable T_ad). The saturation-pressure range check still raises `IndexError`. Results for valid input are bit-identical.

- **`k_haynes230` is the brochure's alpha * cp * rho, not its printed conductivity column (#881).** Haynes brochure H-3000O prints k as a smoothed line (8.9 W/m K at room temperature, +2.0 per 100 C). That line is 13% above the brochure's own diffusivity x specific heat x density at 500 C, and 1.7-2.8 of Georgia Tech Gen3 CSP's 95% measurement bands above its measured k at 340-600 C; the product lies inside the band. `k_haynes230` now interpolates the product, 25-1000 C. k changes by -3% at room temperature, -11% at 500 C (16.33 instead of 18.40), and -6% at 1000 C (26.68 instead of 28.40). Walls of `haynes230` conduct less.
- **Breaking (0.9.0): `k_stainless_steel_316` and the `ss316` material are the Gen3 CSP measured table, valid 533-1093 K (260-820 C) (#881).** The quadratic used before cited no document. It lay inside Gen3 CSP's band at every measured point, but ran from 300 to 1200 K, past the data. Inside 533-1093 K the values move by at most 5%. Outside it, `k_stainless_steel_316(T)` now raises, and `get_material_conductivity("ss316", T)` and a `WallLayer` hold the end value, unless `extrapolate=True`.

- **The CSV export writes the displayed result (#783).** `/solve` keeps its last few solves, keyed by the network's nodes, edges, inputs and solver settings. Layout, labels and the results the GUI stores in each node are left out of the key. `/export` of the same network writes that solve, with the current labels, and does not solve again. It solves only when nothing matches (the canvas changed since the solve, or the backend restarted). The header `X-Combaero-Export-Source` says `displayed` or `solved`, and the GUI shows a notice for a fresh solve.
- **Rows the builder made are labelled (#783).** These are junctions, port chambers, auto links, and an array's rows and segments. Each gets a derived `kind`, a label from its parent's, and a new `parent` column naming the drawn node or edge it belongs to. A drawn row's `parent` is empty. Before, 54 such rows across the palette exported with no kind or label.
- **List cells are strict JSON in the CSV (#783).** `rows_*` and `T_interface` cells are JSON arrays, with a non-finite value as `null`, so `pd.read_csv(path, converters={col: json.loads})` reads them. They used to be Python repr text that wrote `nan` / `inf`. `load_results` reads both forms, including a pre-#783 list holding `nan`, which it used to leave as text.

- **Breaking (0.9.0): outside a material's sourced range, the database follows a k(T) envelope instead of holding the end value (#882).** This applies to `get_material_conductivity(name, T)` and to a `WallLayer` with a named material.
  - The old clamp kinked k(T) at the end of the data, and nothing reported it.
  - Now the end slope continues, ramped linearly to zero over 100 K and then flat. k(T) is C1, and it moves at most |slope| x 50 K past the data, about 1 W/m K for the superalloys.
  - The envelope runs from 200 K, the NASA-9 lower bound (or T_min, if that is lower), up to the material's solidus. Where no solidus is printed, the cap is the generic 1533.15 K (1260 C, the lowest sourced solidus). Alumina, ysz and al6061 are capped at T_max. k is held beyond the envelope.
  - Results change only where a wall layer's temperature leaves its material's range; inside the range they are bit-identical. `extrapolate=True` (linear, unbounded) is unchanged.
- **New flag: `material_in_range(name, T)`, and `"<wall>.material_in_range"` in network results** (0 when any layer's k came from the envelope). The GUI wall panel says so. `get_material_info` also reports `T_envelope_min` / `T_envelope_max`.

- **The C3X example's walls conduct on the log mean of their two faces (#887),** not on the segment's: the metal narrows from the segment to its hole. The gas-side `Nu_multiplier` tuners are refitted (the leading edge goes from 2.2 to 1.9–2.0).
  - Untuned run 158 (cross-source): bias −21.0 → −14.3 K, MAE 26.9 → 23.8 K.
  - Run 149 with run 158's multipliers: MAE 11.7 → 9.6 K, bias +11.2 → +8.9 K.

- **The C3X worked example now runs on burner gas, as the rig did** (Hylton 1983: Allison's natural-gas burner ahead of the cascade). Ambient humid air is compressed to the burner inlet (isentropic efficiency 0.85, assumed). A CH4 burner's fuel flow is iterated so that the single cascade inlet is at the measured Tt1, and the burner's 4 % loss (assumed) leaves Pt1. All 19 segments now draw from that one inlet.
  - With products in place of dry air the tuned multipliers drop by about 2 %. The untuned run 158 error is MAE 26.9 K (bias -21.0 K, worst -83 K at the leading edge), cross-source.
  - Two variants closer to the rig were tried and are documented as left out. One coolant supply temperature for every hole is worse than the measured per-hole temperatures, and points at missing 2-D conduction in the trailing edge (#898). Merging all external paths to one exit pressure is feasible path by path but fragile when solved jointly (#899), and waits for design targets (#897).
  - The burner's fuel flow is a `FuelBoundary` (#892), solved in the network; the script's outer fuel loop is gone, and the GUI example carries the fuel node. The compressor is still computed in the script until #893 exists.

- **`calc_T_from_cp` searches [1, 6000] K, down from [200, 6000] K.** Below the 200 K NASA-9 data the ideal-gas tail (#606) gives cp at every T, so a cp only the tail reaches now returns its T instead of raising `RuntimeError`.
  - **Before:** a target such as H2's cp at 50 K raised "no root in [200, 6000] K", and a `T_guess` below 200 K started the scan at 200 K.
  - **Now:** the root nearest `T_guess` in [1, 6000] K (the message says "no root in [1, 6000] K"). The scan's nodes are 10 K apart, 10 % of T below 100 K, and each step goes to the side whose next node is nearer the guess. The lower end is 1 K because the tail's cp approaches its 0 K value as T^3: at 1 K, T is still fixed to about 5e-7 relative (CO, the flattest species), at 0.1 K only to about 5e-4.
  - **Results that change:** a guess below 200 K, and a target that also has a root in the tail nearer the guess than the one above 200 K (N2's cp from 200 to 340 K is also reached between about 35 and 127 K). For a guess in [200, 6000] K and a target with no root below 200 K the result is bit-identical.

- The project moved to the `thirox-lab` GitHub organisation: the package URLs (homepage, repository, issues, changelog) and the README links now point at `github.com/thirox-lab/combaero`. The old `thiemom/combaero` URLs redirect.

### Removed

- **The 0.8.0 stagnation spellings are gone (#782).** Deprecated in 0.8.0 with a `DeprecationWarning`; an old keyword now raises `TypeError`, an old name or field `AttributeError`.
  - Functions: `T0_from_static`, `T0_from_static_v`, `P0_from_static`, `solve_P0_from_mdot` (Python and C++) and `T0_/P0_from_static_and_jacobian_M` (C++ `solver_interface.h` and `combaero._solver_tools`). Use the `Tt`/`Pt` names.
  - Keywords: `T0=`/`P0=` of `nozzle_flow`, `nozzle_cd`, `nozzle_thrust_cd`, `solve_A_eff_from_mdot`, `solve_P_back_from_mdot`, `solve_Pt_from_mdot`, `critical_pressure_ratio`, `mach_from_pressure_ratio`, `mass_flux_isentropic`, `duct_static_state`, `duct_htc_static_state(_dY)`, `T_from_stagnation`, `P_from_stagnation`; `compressible.nozzle_flow(T=, P=)`; `heat_entropy_flux(T0_in=, T0_out=)`; `combustor_head_loss(T_in=, P_in=)`; `chamber_merge_face_state(T_main=, T=)`; and the internal solver tools' `P_total=`, `T_total=`, `T=`, `T_up=`, `p0=`, `t0=`, `P_total_up=`, `P_total_down=`. Use `Tt=`, `Pt=`, `Tt_in=`, `Pt_in=`, `Tt_up=`, `Pt_up=`, `Pt_down=`.
  - Fields: `MixtureState.T_total`/`.P_total` (and the constructor keywords), `MassStream.T`/`.P_total`, `StreamJacobian.d_T`/`.d_P_total`, `NozzleSolution.T0`/`.P0`, `DuctStaticState.d*_dT0`/`d*_dP0`, `SonicStatic.dP_dT0`, `OrificeResult.d_mdot_dT_up`, `ChannelSolverResult.d_dP_dT_up`, `CombustorLossResult.dxi_dT_in`/`dxi_dP_in`, `MixerResult.T_mix`, `.P_total_mix`, `.dT_mix_*`, `.dP_total_mix_d_stream`. Use the `Tt`/`Pt` spellings.
  - `PressureLossContext.T_in`/`.P_in`: use `Tt_in`/`Pt_in`. A hand-built context without them raises `AttributeError` naming the new fields.
- **Network results: the 0.8 static-T transition aids are gone (#782).**
  - The `"<chamber>.T_static"` result key: use `"<chamber>.T"`.
  - The once-per-process `FutureWarning` on reading a `.T`, `T_in`/`T_out` or `port_<i>_T` key: `NetworkSolver.solve()` returns a plain `dict` again.
  - `NetworkMixtureState(P, Pt, T, Tt, m_dot, Y)` no longer takes `T`, and `state.T = ...` no longer assigns `Tt`. `Tt`, `m_dot` and `Y` are keyword-only, so a positional call written for the 0.8 signature raises `TypeError` instead of shifting its arguments.
  - A custom element's Jacobian keyed `"<node>.T"` (the 0.7 name of the stagnation-T column) raises `KeyError` naming `"<node>.Tt"`; 0.8 relayed it with a warning.
- **The Stolz metering Cd (ISO 5167:1980) is gone (#876).** `MeteringCdCorrelation.Stolz`, C++ `Cd_Stolz` and `OrificeElement(correlation="Stolz")` (now a `ValueError` naming the replacement). Use `ReaderHarrisGallagher`, which ISO adopted in its place. The internal `_metering_cd_and_dRe`/`_metering_cd_in_range` bindings, which only existed to skip the Stolz warning, went with it.
  - GUI: the Inspector's Stolz option and `SolveResponse.deprecations` are gone. A saved file with `"correlation": "Stolz"` is refused (`/solve` returns HTTP 400 naming the orifice and `ReaderHarrisGallagher`), not silently migrated.
- **`k_aluminum_6061` and the `"al6061"` material are gone (#761).** They cited no source. Give the layer its conductivity (`material="custom"`) instead.

### Fixed

- **Cold start of channel chains and compressible elements (#971).** The resistance-network seed now matches what each element's residual does.
  - **Before:** every `ChannelElement` was charged one exit head, which its residual loses only into a static-coupled pressure boundary (#360). The D-7156 film rig's ten 1.22 cm segments started 14 kPa high against a 13 Pa root. All 62 runs stalled cold, at |F| ~ 9, and needed the warm-start retry (10-27 s).
  - **Now:** a channel's seed charges the head only where its residual does. A `CrossflowSegmentElement` always charges it. A compressible channel or orifice whose seed flux reaches 0.3 of the sonic flux has its seed resistance re-fitted to the flow its own C++ kernel (Fanno, or the orifice's compressible model) passes at the seed's pressures. 61 of the 62 D-7156 runs converge cold, in a median 3.6 s. Run 19 (coolant passage at Mach 0.7) still retries.
- **A `FuelBoundary` on an equilibrium combustor starts on the equilibrium curve (#991).** Its start burned the guessed oxidizer by complete combustion whatever the combustor's `method`; rich of stoichiometric that runs hundreds of kelvin hotter (2020 K against 1500 K at phi 2.2), so a rich target started near twice its fuel flow and reached the rich root only through the LM fallback, or landed lean. The start now uses `combustion_equilibrium` for `method="equilibrium"` and lands on the root; complete-combustion combustors are unchanged.
- **Orifice Jacobian: composition columns and the Reynolds-dependent Cd (#964).** Only the Jacobian changes; residuals and solved values are unchanged.
  - **Reynolds-dependent Cd and composition.** Re ~ |m| / mu, and mu depends on the composition. The derivative carried mu's temperature column but not its composition column. This affects Reader-Harris/Gallagher, Stolz, the Idelchik holes, McGreehan-Schotsch and Lichtarowicz. With an RHG orifice downstream of a combustor, the orifice row's fuel-flow column was 7e-4 off past the approach pipe's sonic limit. The Rogers-Hersh crossflow Cd on `ImpingementPlateElement` had the same gap in its density.
  - **Kernel composition (Y) columns.** These were forward differences, up to 3e-3 of the row off a central difference. They are now analytic in all orifice kernels and in `approach_static_pressure`:
    - incompressible kernel: through the density;
    - nozzle kernel: through its isentrope, and choked at a fixed P*;
    - `jobson1955` and `iso5167`: through gamma and R, as two more DualN directions;
    - `approach_static_pressure`: in its low-Mach root, subsonic, and past choke.
  - **Idelchik beveled and rounded holes.** Their tables' geometry slope was returned as dCd/dRe (0.83 and 1.53 per unit Re, for a Cd that does not vary with Re). The self-consistent orifice row then had nearly zero in every column except its own flow. The Re slope is now 0.

- **In-pipe orifice Jacobian past its approach pipe's sonic limit (#958).** Past choke, `duct_static_state` (and `approach_static_pressure`, which an in-pipe orifice drives from) returned the sonic pressure as `critical_pressure_ratio(Tt, Pt) * Pt` with `dP_dTt = 0`. With a variable cp, P*/Pt moves with Tt (18.5 Pa/K for air at 1900 K and 40 bar). The ratio also came from a golden-section search on G(P) that is good to about 1e-8, and central differences read that noise as a derivative. In #788's network, a channel choked into a plenum ahead of an in-pipe exit orifice, the orifice row's air and fuel flow columns (through the combustor's mixing temperature) were 7.2 analytic against 3.4 by central difference. The choked state is now the sonic point on the isentrope through (Tt, Pt): T* from the M = 1 energy balance, P* = Pt exp(-(s(Tt) - s(T*))/R), with analytic Tt and Pt partials, M = 1 exactly, and T, rho and their partials of that state (`dT_dPt` is now 0; it was nonzero). P* moves by about 4e-8 relative. The subsonic branch and the choke decision are unchanged.

- **`jobson_expansibility` converges at high beta Ci (#957).** It threw at its iteration cap on valid input: beta 0.5-0.9 with Ci 0.6-1, over a band of pressure ratios, 2747 of a 33000-point sweep. All of those points are outside Jobson's range (Ci <= 0.7, beta <= 0.75). The approach-velocity consistency was a fixed point on the jet flux. It cycles wherever the flux map's slope at the root is below -1. Now it is a bracketed Newton (rtsafe) on the same equation, with the flux slope by forward mode. The root is unique. Every sweep point converges.
  - **The former failures** agree with a 40-digit reference to 1e-13 for r < 0.9. Near r = 1 (r >= 0.966) they agree only to 3.8e-12: the flux itself carries that rounding there.
  - **Where the old iteration converged**, values move by up to 2.3e-13 relative for r < 0.9, and by up to 3.7e-11 near r = 1. There the new value is the closer one to the reference in 77 of 81 cases.
  - **Not changed:** the network orifice kernel never called this function.

- **Static temperatures below 10 K are solved, not floored (#935).** The isentropic and stagnation static-T solves held their iterate at `T_FLOOR` = 10 K and, when the root lay below it, returned 10 K, which is not a root. The mass flux and Mach built on that value were the floor's. This happens for an expansion to P/Pt below ~1e-5 from 300 K, for `T_from_stagnation` at hypersonic M, and for `Tt_from_static` from a static T below 10 K. The gas is ideal to 0 K (#626), so these solves now return the actual root, however small. From 300 K air, P/Pt = 1e-8 now gives about 1.5 K and the matching flux and Mach. The affected functions are `mass_flux_isentropic`, `mach_from_pressure_ratio`, `T_from_stagnation`, `P_from_stagnation`, `Tt_from_static`, `Tt_from_static_v`, `Pt_from_static`, `stagnation_from_static` and the Fanno-march energy solve. Below 10 K the isentropic Newton steps in ln T, so it reaches a root many decades down within its budget. Values whose solve never reached the floor are bit-identical.
- **Entropy is exact below 10 K (#935).** `s`, `s_mass` and `s_R` held their value below 10 K, so s was flat there while ds/dT = cp/T was not. They now follow cp ln T down to 0 K and hold only below 1e-300 K, where a non-positive solver probe stays finite. The low tail also returned NaN below ~1e-14 K, and it no longer does. `h_RT` and `g_over_RT`, which feed the equilibrium constants, still hold below 10 K.

- **The Jacobian of a `ThermalWall` with k(T) layers now includes the conductivity's temperature dependence (#932).** A layer of a database material takes its k at its mean temperature in the wall profile, iterated to a fixed point. The Jacobian held that k fixed, so dQ/dh and dQ/dT_aw missed the part of the heat that moves with the wall temperatures through k. The worst scaled Jacobian error on a two-duct network was 3.5e-5 to 3.5e-4 against central differences for ss310, a YSZ + bond coat + Inconel 718 stack and Hastelloy X past its data. It is now 2e-8 to 5e-7, the same as a constant-k wall under the same conditions. The fixed point is differentiated implicitly in C++, and dk/dT comes analytically from the material data, including the #882 envelope's fade. Residual values are unchanged bit for bit; only the Jacobian of k(T) walls changes.
- **Added: `get_material_conductivity_and_slope(name, T, extrapolate=False)`.** It returns `(k, dk/dT)`, with the same k as `get_material_conductivity` bit for bit. The slope is a table segment's slope inside the data, the envelope fade's slope outside it, the end slope with `extrapolate=True`, and 0 where k is held or floored. In C++ it is `materials::material_conductivity_and_slope`.
- **Added: `wall_resistance_fixed_point_derivatives` and `WallCouplingResult.dQ_dR_wall`.** The first gives the total derivatives of a k(T) wall's conductive resistance at its fixed point. The second is dQ/d(t/k) of any one layer.

- **Equilibrium inner solves reach their roots (#927).** `wgs_equilibrium`, `reforming_equilibrium` and their adiabatic forms, and so `combustion_equilibrium` and the equilibrium `CombustorNode`, used a damped Newton with a finite-difference slope that stopped at its iteration cap without saying so.
  - **Before:** with CO or H2O a trace, the 1e-6 step was wider than the admissible WGS interval and the extent stopped at ln(Q/K) = 12-14. Where the equilibrium lies away from the trace, a whole shift was missed: CO came back at 7e-7 instead of 0.015 for CO 1e-8 at 800 K. The adiabatic temperature cycled on a rich reforming case (a burned CH4 state at 1500 K) and returned 1390 K, reported converged, against a root at 924 K (Cantera HP on the same species: 923.85 K). An unreacting inlet colder than 300 K came back at 300 K, the old bracket end.
  - **Now:** the WGS extent is a bracketed Newton with the analytic slope, converged to |ln(Q/K)| <= 1e-12 or to its rounding error. The adiabatic temperature is bracketed from its guess and solved by Brent's method on the energy residual, between 10 K and 4500 K. An empty interval returns 0 straight away. An iteration cap, a NaN, or no root in the range now raises instead of returning the last iterate. `reforming_equilibrium_adiabatic` reports `converged = False` when an inner equilibrium at either end of the final temperature bracket did not converge; before, it dropped that flag.
  - The suite's 787 cap exits (and 42 in the Cantera validation) are now 0. Every other result in the suite and in the Cantera validation is unchanged to 5e-12 K and 3e-15 in mole fraction.

- **A ribbed or pin-fin channel fed straight from a reservoir takes its pressure drop at the duct's static density (#926).** Ribbed and pin-fin surfaces own a `ChannelElement`'s drop and read the from node's (Tt, P) for it. Where that node holds its gas at rest (`state_at_rest`: `PlenumNode`, `CombustorNode`, `WallNode`, `PressureBoundary`, `MassFlowBoundary`) its state is the reservoir's, P = Pt, so the density was rho(Tt, Pt) and not the duct's.
  - **Before:** the density was high by about M^2 / 2 (0.5% at Mach 0.1, 4.6% at 0.3, 13% at 0.5, 26% at 0.7, 46% at 0.9), and the drop at a given flow low by 0.5%, 4.4%, 11.5%, 21% and 31%. A ribbed duct from a 2.6 bar reservoir to 2 bar passed 6.3% too much air (0.0379 against 0.0357 kg/s, inlet Mach 0.50).
  - **Now:** at such an end the drop reads the duct's own inlet static state, the isentropic acceleration from the node's (Tt, Pt) to the duct's mass flux (`duct_htc_static_state`, the state its heat transfer already used, #636), for the density and for the viscosity in the pins' Reynolds number, with analytic partials in the flow and the node's Tt, Pt and composition. `mach_in` and the inlet reference block report that state too.
  - **Unchanged:** every other channel and network (F and J hash-identical over the 520 stress networks without a reservoir-fed pin-fin channel), ribbed and pin-fin channels leaving a momentum chamber (Pt > P), an outflow into a boundary that takes it at its static pressure, `CrossflowSegmentElement`, and the cooling scorecards (they score the correlations directly).

- **GUI: a flow edge that names a handle its node does not have is reconnected on load, or reported (#901).** The canvas draws an edge only between handles its two nodes define, while the backend solves the edge by its endpoints alone. A file with a wrong handle id, typically one written by a script, therefore solved but showed the node disconnected.
  - **Interchangeable handles** (a plenum's four sides, an unnamed handle): the edge end moves to the node's first handle on that side.
  - **Named ports** (an ejector's secondary inlet, a chamber's side inlet): the end is not guessed. The load notice lists it to reconnect by hand.
  - **No handle id at all** is valid and left as it is: the canvas uses the node's first handle on that side.
  - **Handle table:** the ids come from a generated table of every node type's handles, `nodeHandles.generated.json`, written by `scripts/gui_node_handles.py` and kept current by a test.

- **The cold start no longer mixes unrelated parts of a network (#899).** The reference state that seeds a cold solve (pressure, temperature, composition, mass flow) is now computed per flow component: the parts no element connects. Before, one global reference averaged every boundary. A disconnected 358 K coolant duct beside 19 near-sonic 808 K gas paths pulled the gas reference to 583 K, and every gas path started 18 % high. The solve converged, but with nine momentum chambers on the supersonic root of their closure, so it was reported inconsistent. It now converges on the subsonic roots, with the flows the gas paths have without the duct.
  - **Temperature is shared within a thermal group.** The reference temperature is taken per thermal group: the flow components joined by thermal walls. A stream heated through a wall ends up between its own boundary temperatures and its partner's. Taken from its own boundaries alone, the cold stream of a hot/cold pair started 60 % too dense and solved 4-7 times slower.
  - **A network of one flow component starts exactly as before.**

- **A `ChannelElement`'s end state at a node at rest is the duct's own (#888).** Where a channel ends on a node whose own definition holds Pt = P (new class attribute `state_at_rest`: `PlenumNode`, `CombustorNode`, `WallNode`, `PressureBoundary`, `MassFlowBoundary`), `mach_in` / `mach_out` and the inlet reference block (`velocity`, `Re`, `rho_ref`, `mu_ref`, `Pr_ref`, ...) used to be read at the node's state: a reservoir (P = Pt, T = Tt) where the duct leaves it, or a static pressure with the arriving Tt where the duct enters it.
  - **Before:** a 1 mm duct from a reservoir (Tt 808 K, Pt 243 kPa) to P/Pt 0.6 reported `mach_in` 0.577 and `mach_out` 0.962, although the duct carries Mach 0.895 end to end. At P/Pt 0.537 the unchoked duct reported `mach_out` 1.085. In the C3X vane example the gas segments, which leave the cascade-inlet plenum, reported `mach_in` 0.57-0.58 on the suction side, and eight of them reported `mach_out` > 1.
  - **Now:** a compressible duct reports the Fanno station at that end's Mach, so `mach_in` / `mach_out` equal `M_in_duct` / `M_exit_duct` there: 0.895 / 0.896 in the 1 mm duct, and 0.86-0.95 at the C3X suction-side inlets. An incompressible duct reports the residual's own density state (Tt at the approach static P) where the node holds a total pressure. The second-law check reads these Mach numbers, so its `sigma_norm` moves at such ends.
  - **Unchanged:** every residual and Jacobian (hash-identical over 560 stress networks), every solution, `T_aw` and `htc` (the correlation was already at the duct's inlet static state, #636), ends on a momentum chamber (Pt > P), `P_in` ... `Tt_out`, ribbed and pin-fin surfaces, and `CrossflowSegmentElement`.

- An `OrificeElement` in a pipe (`in_pipe=True`, a velocity of approach) no longer stops a solve on a spurious reversed flow (#875). Its row solves its own flow consistently with the approach pipe's static pressure. When an outer iterate overshot the flux at which that static pressure falls to the outlet's, that inner solve could run down the mirrored branch to the approach pipe's choke. There it stopped at a reversed flow tens of times the true one, where the row has no pressure slope. It is reachable only when C E beta^2 > 1, i.e. the jet's effective area exceeds the pipe's: a fixed Cd 0.98 at beta 0.95, but also Reader-Harris/Gallagher's Cd near 0.6 at beta 0.95. The inner solve now keeps the root on the side of its drive. Rows that already reached that root are bit-identical.
- `combustor_liner_cooling_network.py`: its two exits (fixed Cd 0.98, bore 0.95 D) are back in their pipes (`in_pipe=True`, as the example was written), undoing 0.8.0's `in_pipe=False` workaround now that #875 is fixed. Its printed results are unchanged. In the coupled case the coolant exit orifice's Pt drop falls from 14.88 to 14.48 kPa, the velocity of approach (E) now applying.

- **Inner solves report a cap exit instead of returning a non-root (#868).** An audit of every iterative solve in the C++ kernels (52 loops; table in `docs/archive/AUDIT_868_INNER_SOLVES.md`) found loops that, on reaching their iteration cap, returned their last iterate as if it were a root.
  - **Before:** `calc_T_from_cp` returned the 6000 K clamp for a target cp no T in [200, 6000] K has; `fanno_inlet_mach` was up to 73 % high below ~1e-10 of the sonic flux (its slope stencil crossed M = 0 and Newton 2-cycled); the Fanno/Rayleigh choke search ran Brent on an interval with no sign change and returned a point of it as the choked flux (a probe at an unphysical exit Tt in the C3X vane example); `chamber_static_state` returned T = 1e11 K for a probe at Tt = -1.4e19 K; a caller's small `max_iter` in `friction_colebrook`, `Tt_from_static` or `T_from_stagnation` returned the half-converged iterate.
  - **Now:** these raise `RuntimeError` (C++ `std::runtime_error`) naming the function and its last residual; `fanno_inlet_mach` returns the root. Bracketed solves still return a cap exit whose bracket or last step has closed to 1e-10 (a root to rounding). The 10 K `T_FLOOR` clamp of the static temperature is unchanged.
  - Results are bit-identical where the solves converge: the network F and J at the start and three perturbed points of all 14 stress families x 40 seeds (2240 evaluations) and sweeps of the touched kernels hash identically to main. A network whose solver probes reached one of these exits rejects that probe as an unphysical step instead; the C3X example's tuned run lands on the same root within 1e-11.

- **A momentum chamber's closure has one root per total pressure (#788).** At a fixed flux the compressible closure Pt(P) has its minimum Pt* at the sonic static pressure, so one Pt had a subsonic and a supersonic static root. Newton could land on the supersonic one, as in the issue's channel-to-orifice junction at Mach 1.46. Below the sonic point the closure is now continued by its mirror image about Pt*. That continuation is monotone and C1 through the sonic point, with analytic partials. Seeded on the old supersonic root, the issue's 95 mm network now finds its subsonic root (Mach 0.74). Before, it ended INCONSISTENT at Mach 1.15. Subsonic evaluations are unchanged bit for bit. Over the 14 stress families x 40 seeds, F and J are identical at x0 and three perturbations in 558 of 560 networks. The two that differ are the two whose chambers are evaluated past the sonic point there.
- **A choked channel discharging into a junction now has a root: the choked-bounded junction (#788).** A choked duct into a pressure boundary had a root, because the boundary takes the flow "choked bounded" (#636). Into a junction it had none: no_progress at exit bores 97-99.5 mm in the issue's network, and with a bore equal to the pipe or a lumped orifice wider than the pipe. When the downstream element needs a Pt below the sonic Pt* of the pipe's flux, the junction's root lies on the continuation. The solver accepts that root as `__choked_bounded_junctions__` when one choked compressible channel of no larger area feeds the junction and no other row reads the junction's static pressure. That root is the one a `PlenumNode` in the junction's place gives. The combustor's pressure equals the channel-straight-into-the-outlet network's to 1e-7. The expansion loss Pt* - Pt sits in the channel's control volume, where the second-law check passes, and mass and energy close. The junction reports as that plenum: `P` = `Pt`, `mach` 0, plus `choked_bounded`, `P_sonic`, `Pt_sonic` and `expansion_loss`. The channel reports its choked exit (`mach_out` 1).
- **Other roots past a chamber's sonic limit are refused with the reason (#788).** This covers a chamber fed by anything but one choked channel, a channel wider than the chamber, and an element reading the chamber's static pressure. These roots stay in `__supersonic_chambers__` and are reported INCONSISTENT, now with the reason in `__message__`. Two stress draws that stalled (ejector 32, random 26) now end on such a refused root.
- **A chamber that starts past its sonic point is seeded on its closure's root (#788).** The closure's slope is zero at the sonic point. The start P = Pt of a bore-equal-to-pipe junction lay past that point, and the first Newton step crossed the fold to a negative P. That cold solve stalled without the Newton-homotopy retry; it now converges cold. No other start moves: F and J at x0 change only in the two stress networks whose start is past the sonic point.

- **More NaN input refused instead of answered (#769 follow-up).** Replaying every numeric call the test suite makes with NaN in each argument found arguments where a clamp, a root finder or a "no correction" branch turned NaN into a plausible finite number. They now raise `ValueError`; finite values keep their existing meaning.
  - **Fanno and Fanno-Rayleigh:** `roughness` and `f_multiplier` in `fanno_duct`, `fanno_choked_mass_flux`, `fanno_length_between`, `fanno_dx_dmach`, `fanno_channel_flow`, `fanno_rayleigh_*` (a NaN roughness or f_multiplier still returned a finite choked flux); `Pt0`, `Tt0`, `Tt_exit` in `fanno_rayleigh_choked_mass_flux`; L and D must be finite.
  - **Orifice:** `velocity_of_approach_factor(nan)` returned 1; `expansibility_factor` kappa; `jobson_force_defect` / `jobson_ci_from_cd` beta; `Cd`, `Cv`, `area` and `beta` of the Jobson and ISO 5167 orifice residual kernels.
  - **Others:** Baldauf film effectiveness `Tu` (the `Tu <= 0` check let NaN through); `T_adiabatic_wall` T_static and v (the 10 K floor answered 10 K for a NaN T_static at v = 0); area-change F_in, F_out, D_h and the conical length; `impingement_eta_r` xn_d and yn_d (they set the extrapolation flag); `omega22_and_derivative` delta_star (clamped to 0); `combustor_fraction_loss` m_ref.
- **Ejector kernels validate their input (#769).** `src/ejector.cpp` checked nothing: a NaN area ratio gave a finite residual and `ejector_huang_mixing_coefficient(nan)` returned 0.84. Every public ejector function now raises `ValueError` for non-finite input and for parameters outside their physical domain: gamma > 1, r_gas > 0, nozzle efficiencies in (0, 1], recovery in (0, 1], 1 < A_p1/A_t <= A_3/A_t, positive areas, A_exit >= A_throat, loss coefficients and eps_frac >= 0, blend windows lo < hi. The standalone closures (P_c*, P_b0, entrainment, choked and C-D nozzle flow, jet-pump discharge) also need a finite, positive supply stagnation state. `ejector_element_residuals_and_jacobian` and `ejector_jetpump_incompressible_residuals_and_jacobian` check their Newton iterates for finiteness only, since a trial point can carry a negative pressure or a reversed flow. `EjectorElement` refuses NaN and inf areas, recovery and loss coefficients at construction; its checks were written `x <= 0` and let NaN through.
- **No neutral answer from NaN (#769).** NaN no longer lands in a "no model" or "no flow" branch and answers its neutral value. That branch is the "no fin model" eta = 1 in `pin_fin_array_efficiency` (+inf `k_pin` and +inf h stay accepted as limits), Gc/Gj = 0 in `crossflow_to_jet_ratio_at_x`/`_at_row`, or Y = 1 (or 0) in `mcgreehan_schotsch_1988_expansion_*`. It is also zero flow in `fanno_channel_flow`, `fanno_rayleigh_channel_flow`, `area_change_loss`, `jet_impulse`, `crossflow_velocity_ratio`, `rogers_hersh_1976_cd`, `florschuetz_metzger_su_1982_cd`, `combustor_fraction_loss` and the Jobson and ISO 5167 orifice residuals. These functions now raise `ValueError`; finite input, including the zero, negative and gamma <= 1 conventions, answers as before.
- **Ejector `recovery_efficiency` is capped at 1 (#769).** It multiplies the lossless mixed stagnation pressure, so a value above 1 is a pressure gain beyond lossless mixing. It is now refused with a `ValueError` naming it, by `EjectorElement`, the C++ ejector closures and both element kernels. The GUI enforces the same bound: the `EjectorData` schema field (`0 < recovery_efficiency <= 1`), the frontend validation and the Inspector input. No shipped example, GUI example or test used a value above 1. Previously any positive value was accepted.

- **`chamber_static_state` refuses a non-positive or non-finite stagnation T (#937).** A solver probe whose upstream state mixes a garbage chamber Tt (-4e18 K in a reversed branching-tee manifold) reached the static-state Newton.
  - **Before:** Tt = -1 K returned a static T of -1.00001 K and Tt = 0 returned 0; larger garbage ran the Newton to its cap.
  - **Now:** `ValueError` ("chamber_static_state: Tt must be finite and > 0 K"), an invalid argument as in #928, which the network solver rejects as an unphysical probe. Physical states are unchanged.
- **The ejector's two bisections flag a target outside their bracket (#945).** Flag and clamp: the value is still the bracket end, now reported.
  - **Before:** the nozzle-exit Mach (area ratio outside A/A* of [1, 50]) and the C-D nozzle's exit-static inverse (m_dot above the flow at P* or below the flow at Pt (1 - 1e-12)) returned the bracket end silently, the inverse with implicit-function partials of a root that does not exist there.
  - **Now:** `EjectorCDExitStaticJacobian.out_of_range` (+1 above the choked end, -1 below the near-zero end) with the clamped function's partials (0, k, 0); the new internal `ejector_nozzle_exit_mach` and `EjectorElementResidualJacobian.nozzle_exit_mach_out_of_range` carry the Mach flag, and `EjectorElement` reports it as `<id>.nozzle_exit_mach_out_of_range` (0 in range; the GUI result lists it). Neither bracket is left at any EjectorElement iterate of the stress families; in-range results are bit-identical.
- **Cap rule for `duct_static_state` (#941).** Its bracketed ln P Newton now applies the #868 cap rule (a cap exit whose last step is not closed to 1e-10 raises `RuntimeError`) instead of returning the last iterate. Never reached; results bit-identical.
- **Removed the dead `solve_station_state` (#942)** from `src/compressible.cpp` (an internal, `[[maybe_unused]]` helper with no caller).

- **Five inner solves reach their roots exactly (#936, #938, #939, #940, #943).** Follow-ups to the #868 audit; each moves results slightly.
  - **`T_from_stagnation` / `P_from_stagnation` (#936):** the Newton slope keeps the T dgamma/dT term, so it converges quadratically. The energy residual falls from up to 2.4e-8 cp T to 9e-14 cp T; static T moves by up to 1.1e-9 relative and static P by up to 5.7e-9.
  - **`State.set_VH_mass` (#938):** it is the enthalpy inverse `calc_T_from_h_mass` (u + R T = h), as `set_HP_mass` is, in place of a hand Newton that stopped at 1e-4 K. T moves by up to 2.2e-8 K over the sweep.
  - **`calc_T_from_cp` (#939):** returns the root nearest `T_guess` in [200, 6000] K. cp(T) is not monotone (N2 has a minimum near 214 K), so the guess picks the branch. The solve scans outward from the guess, clamped into the range, in 10 K steps, then runs a bracketed Newton. The range is the span where every species has NASA-9 data, because the tails level cp off beyond it. No root in the range raises `RuntimeError` ("no root in [200, 6000] K"). A non-finite or non-positive `cp_target`, `T_guess` or `tol` raises `ValueError`; before, a negative guess returned 530 K for air. Roots move by up to 8e-7 K.
  - **`fanno_inlet_mach` (#940):** Newton uses the exact flux slope in place of a finite-difference stencil. M_in moves by up to 3.4e-14 relative.
  - **`jobson_expansibility` (#943):** the fixed-point tolerance is 1e-13 (`inner_solve::kJobsonExpansibilityRelTol`), not 1e-15. The old tolerance cycled in the last bits. Results move by up to 9.8e-14 relative.

- **A fuel controller's flow no longer sets the cold-start reference flow (#899).** A `FuelBoundary` is a mass-flow boundary whose flow is solved, not given, but its `m_dot` (whatever the last start left in it) was added to the reference flow of its flow component. Without it, the C3X burner network's cold GUI solve takes 15 s instead of 21 s, the same as with a fixed fuel flow. Networks without a fuel controller start exactly as before.
  - The issue's broader proposal was tried and dropped: estimating the reference flow from the pressures wherever a component has two or more pressure boundaries, so that mass-flow boundaries seed only their own paths. On the stress networks it lost one random seed (27) and gained none.

### Security

- **GUI frontend lockfile refreshed against `pnpm audit`.** The audit found 12 advisories (7 high), all in transitive packages that Dependabot's grouped updates do not move; 1 moderate remains.
  - Shipped-dependency side, now clean under `pnpm audit --prod`: follow-redirects 1.16.0 → 1.16.2 (under axios) and source-map-js 1.2.1 → 1.2.2.
  - Dev tooling: zod 4.3.6 → 4.6.7, brace-expansion 5.0.9 → 5.0.12, @babel/core 7.29.0 → 7.29.7, and vitest ^4.1.11 → ^4.1.12 (which brings @vitest/mocker 4.1.12).
  - Left as is: minimatch 10.2.6 under eslint (a moderate ReDoS). Even the newest eslint requires minimatch ^10, and the fix is in 11. It only globs this repository's own files.

### Documentation

- **`dewpoint`'s RH range documented (#944).** RH in (0, 1]; RH = 1 returns the dry-bulb T; a supersaturated RH > 1, RH <= 0 or a non-finite RH raises `ValueError` (the vapour-pressure step refuses RH > 1 before the dewpoint solve, which is clamped to T and would otherwise return T). In the C++ header, the binding's docstring, `docs/API_CPP.md` and `docs/API_PYTHON.md`; no behaviour change.

- **`k_inconel625` records that Gen3 CSP measures 14-25% more (#881).** The gap is at 260-700 C and comes from the diffusivity. The Special Metals table (Battelle, annealed) is kept until a third source settles which is right.

- **What a `ThermalWall` represents is documented** (class docstring and `docs/API_PYTHON.md`). It is one shared patch through which heat crosses normal to the wall. Its area is the shared part of each perimeter. An element may carry one wall per shared face, and faces without one are adiabatic. Each side uses its element's perimeter-mean h. No heat conducts along the metal between elements in series (#898).

## [0.8.0] - 2026-10-31

### Highlights

The 0.8 release of the breaking changes announced in 0.7 (61 commits since
0.7.1); the sections below are the record. What changes numbers or code on
upgrade:

- **`T` is static, `Tt` total, everywhere (#782).** Network results report the
  static temperature under `T`, as they already did for `P`; stagnation inputs
  are spelled `Tt`/`Pt` across the Python and C++ API. The old spellings warn
  and are removed in 0.9.0. Four kernels that needed a static T now get one.
- **Results tables change shape (#774, #784).** Every CSV/`to_dataframe()`
  quantity column carries its unit (`<key> [<unit>]`), exports carry the
  case's inputs, and `NetworkRunner.sweep()` returns one row per case with
  `(kind, entity, quantity)` columns. `load_results()` reads old and new.
- **Orifices (#788, #861, #872, #873).** `jobson1955` is the default
  compressible model for incompressible-referenced Cd; a bore against its pipe
  follows a rule instead of a silent wall-hole fallback; `Cd` is accepted only
  with the `"fixed"` correlation; Reader-Harris/Gallagher is held below its
  Re_D bound. Stolz is deprecated (#876).
- **Impingement and crossflow (#800, #802, #813).** Jet Cd depends on the
  crossflow by default, crossflow stations take the jets' axial momentum, and
  segments run Fanno from the true stagnation state.
- **The 0.7 deprecations are removed**, among them `BorderCarnotLossElement`,
  `VortexElement(r_c=, n=, omega_rpm=)`, `smr_wgs_equilibrium` and
  `init_strategy='incompressible_warmstart'`.

New: worked examples from papers (scripts and GUI networks), a set-outlet-T
heat source, tapered crossflow segments, opt-in compressible orifice models,
sourced liner/transition alloys with opt-in conductivity extrapolation, and a
stall handover for far-from-root crawls.

### Added

- **Opt-in extrapolation for material conductivity** (#761). The table-backed `k_*` functions, `get_material_conductivity` and `WallLayer` take `extrapolate=False`.
  - **Default (unchanged):** outside its valid range a `k_*` function raises, and the database lookup holds the end value.
  - **`extrapolate=True`:** continues the end slope linearly. On a table this is the end segment; on a fitted polynomial it is the tangent. It raises if that reaches k <= 0 or T <= 0.
  - The continuation is unsourced; the user asks for it explicitly.
  - `k_tbc_ysz` keeps its own range.

- **Liner and transition-piece alloys, and a TBC bond coat, in the materials database.** Each is a function and a `get_material_conductivity` / `list_materials` entry.
  - `k_hastelloy_x`, `k_haynes188` and `k_haynes263` (300-1273 K) interpolate the manufacturers' tables (Haynes International) linearly, so they are exact at every printed point.
  - `k_bondcoat_nicocraly` (848-1373 K) is the NiCoCrAlY line from Zhu & Miller 1999 (NASA/TM-1999-209069, Fig. 2), limited to the range it is drawn over.
  - `k_haynes25` (298-1273 K, Haynes), `k_inconel625` (116-1255 K, Special Metals) and `k_nimonic75` (293-1273 K, Virgamet; its 20 C value matches Special Metals) interpolate the published tables the same way.
  - `k_copper` (high purity, 100-1200 K) and `k_alumina` (dense Al2O3, 100-2000 K) are the NSRDS-NBS 8 recommended values (Powell, Ho & Liley 1966).
  - `k_x22crmov12_1` (1.4923, 293-923 K) is the straight line between the data sheet's 20 C and 650 C values.
  - These are convenience data. A user with qualified properties supplies their own conductivity.

- **Results table for solved networks, one row per case (`gui.backend.results`, #784).** `load_results(path_or_df)` reads a GUI CSV export, a `to_dataframe()` frame, a `NetworkResult`, a table saved with `df.to_csv()` or a 0.7 full-detail sweep frame into one shape: columns a 3-level MultiIndex `(kind, entity, quantity)` with kind `inputs` or `outputs`, so `df["outputs"]["liner", "Tt"]` is a quantity across cases. `append_case(df, result, name=None)` adds a case and reads the full network's inputs from the result itself; `varied_inputs(df)` returns the inputs that differ between rows.
  - **Inputs** are what defines the case: boundary conditions, geometry, element and mixture settings, and the network-wide physics settings under `_network`. Initial guesses, the warm start and the solver options are not inputs; they are kept in `df.attrs`.
  - **Units** travel in `df.attrs["units"]`, keyed by column, so quantity names stay unit-free. The solve status is under the `_solve` entity of `outputs`.
  - Entities are named by label, falling back to the id; two entities with the same name are refused with a message naming both ids. Inputs are the values the network file gives (a field left at its default is absent); an unlabelled entity is named by its id; `_solve` and `_network` are reserved names.
  - **Every numeric input carries a unit** (`INPUT_UNIT_MAP` / `input_unit` in `gui/backend/units.py`: lengths `m`, areas `m^2`, angles `deg`, coefficients and ratios `-`, ...); text inputs are bare. A test fails on any numeric schema field without one.

- **A heat source with a set outlet temperature, and a GUI Heat Source node (#792).**
  - **`EnergyBoundary(Tt_target=...)`** sets the TOTAL temperature of the plenum, momentum chamber or combustor it is on; the heat flow becomes a solver unknown, `{id}.Q` [W], closed by one temperature row. Its Jacobian column is the existing C++ `dT_mix/dQ`, relayed downstream. On a combustor it is the heat taken out (or put in) after the reaction to reach that exit Tt; the fuel flow is unchanged.
  - **Refused (`ValueError`):** `Tt_target` together with a nonzero `Q` or `fraction` on one boundary, two `Tt_target` on one node, and a `Tt_target` on any other node type. A node with no through-flow does not succeed: its heat is undetermined, reported as NaN and listed in `__heat_undetermined__`.
  - **Reported in every mode:** `{id}.Q` per boundary (set, solved, or its share of the fraction heat) and `{node}.Q`, their sum. The second-law check and the balance ledger now count the solved heat.
  - **GUI:** a Heat Source node, wired by one edge to a plenum, momentum chamber or combustor, sets Q (default 1 kW), a fraction, or an outlet Tt (total). The Inspector shows the heat applied (the solved Q in the outlet-Tt mode) beside the node's total temperature.

- **Worked examples from papers, as scripts and as GUI networks** (#793). Each rebuilds a published case, solves it and compares with the source's own values, labelled per the validation policy.
  - **Florschuetz, Truman & Metzger 1981 jet array** (`python/examples/papers/florschuetz1981_jet_array.py`): the ten-row crossflow build-up against the paper's Eq. 8 and row-by-row Nu/Nu1 against the measured Fig. 6 points. Fidelity: 2.1% mean error at B(5,4)I z/d 2; +4.0% bias and 6.9% mean error over all 242 points (`--all`). The script reads the committed digitised data, so it runs from a repo checkout.
  - **ISO 5167-2 orifice metering run** (`python/examples/papers/iso5167_orifice_metering.py`): the flange-tap differential against the permanent loss at beta 0.3/0.5/0.7, checked against the standard's own equations (C, flow equation, 5.4.1 loss). This is fidelity to the equations; there are no measurements.
  - **GUI:** a new **Examples** menu loads the same networks (`gui/examples/*.json`, shipped as combaero-gui package data, served by `/examples`). Each opens with a dismissible panel giving the source, what is reproduced, where to look and the validation label.

- **Tapered crossflow segments (#812).** `CrossflowSegmentElement(area=A_from, area_to=A_to)` models a tapering annulus or a growing liner. Each end's station half-drop, joining jets and static-to-total closure use that end's own section, and friction uses the mean section. Total pressure carries across the segment apart from friction, so the area change is isentropic.
  - **Incompressible:** a from-end station adds Bernoulli's area term, C++ `segment_area_change_drop`, with analytic partials.
  - **Unchanged:** without `area_to`, or with `area_to` equal to `area`, a segment is bit-identical.
  - **Checks:**
    - A frictionless taper passes Bernoulli's flow to 3e-3.
    - The two regimes agree to O(M^2) at low Mach.
    - Tapered chains pass the second law in both directions. A constant section standing in for the taper used to read sigma_norm of about -0.06.

- **Two opt-in compressible orifice models** (#854). `OrificeElement(compressible_model=...)` picks how the compressible regime expands the flow; the GUI Inspector has the same selector. The default (`"nozzle"`, the isentropic nozzle) is unchanged.
  - **`"jobson1955"`:** Jobson (1955), Proc. IMechE 169, 767. A sharp-edged jet widens as the pressure ratio falls, so the flow keeps rising past choke instead of going flat. Its only input is the incompressible contraction Ci from the element's Cd correlation, with `Cd = Cv * Ci` (`Cv` defaults to 1). In a pipe, it adds an approach-velocity term that is not in the paper: f = 1/Ci - 1/(2 Ci^2) - Cv beta^2 + Cv^2 beta^4/2. It works subcritical and choked, and is C1 through the critical ratio. It reproduces Jobson's Table 2 to the printed digits. Cross-source, it lies inside ISO 5167-2's stated expansibility uncertainty for beta 0.2-0.75 and p2/p1 0.75-0.98.
  - **`"iso5167"`:** ISO 5167-2:2003 Eq. (1) with Eq. (5)'s expansibility. It needs the Reader-Harris/Gallagher correlation. Below p2/p1 = 0.75, Eq. (5) is extrapolated and the `model_in_range` diagnostic reads 0.
  - **New functions:** `jobson_force_defect`, `jobson_ci_from_cd`, `jobson_contraction_coefficient` and `jobson_expansibility`. The solver kernels are `orifice_jobson_residuals_and_jacobian` and `orifice_iso5167_residuals_and_jacobian`, which return an `OrificeModelResult` with analytic derivatives.
  - **New diagnostics:** `compressible_model`, `Ci`, `contraction_C`, `expansibility`, `model_in_range` and `Cv`.

- Validation: Schiller's (1933) superheated-steam orifice tests, digitised from Jobson (1955) Fig. 5b, score `jobson1955` as a second choked series (cross-source, 9 points: Jobson RMS 1.9 % vs the isentropic nozzle's 19.8 % at the same Cd) (#861).

- **Stall handover can be turned off (#864).** `NetworkSolver.solve(stall_handover=False)` turns off the stall handover for that call. Without it, a cold rung does not cut itself short when it stalls: hybr does not hand over to the LM fallback early, and the LM fallback does not hand over to the warm-start retry early. Each cold rung runs to its own share of the budget. The retries and nested solves of the call are covered too. The default (`True`) is unchanged.
  - **GUI:** a "Stall handover" checkbox in Solver Configuration, saved with the drawing. A file saved without the setting loads with it on.
  - **Results table:** the setting is stored with the solver options as `solver:stall_handover` (#784). It is not an input.

### Changed

- **`k_tbc_ysz` (plasma-sprayed) follows Zhu & Miller 2000, NASA/TM-2000-210237 (#761).**
  - **As sprayed:** k = 1.0 W/m K at every temperature, the paper's "about 1 W/m-K". It used to rise from 0.81 to 1.39 W/m K, which the literature does not support.
  - **Sintering:** their Eq. 1, with t in seconds. At 1320 C the conductivity rises about 44% within hours, against +3.8% after 1000 h before. There is no 1073 K on/off gate any more, so k is continuous in T.
  - **Fully sintered limit:** k_inf = 1.74 W/m K, inferred from their Fig. 4 plateaus and flagged as unconfirmed.
  - **Caveat:** the measured k is an effective through-thickness value for a 254 um coating under a strong gradient. Users with a characterised coating should supply its conductivity.
  - **EB-PVD:** keeps its as-deposited fit, now marked uncited, with no sintering.

- **Breaking (0.8.0): every quantity column of the CSV export and of `NetworkRunner.to_dataframe()` / full-detail `sweep()` now carries its unit (#774).** About 130 result keys used to have a unit in the GUI only and kept a bare column; they are now labelled `<key> [<unit>]`, which renames them. Code that reads one of them by name must use the new header, e.g. `df["entrainment_ratio"]` becomes `df["entrainment_ratio [-]"]`.
  - **Renamed:** every key formerly in `GUI_ONLY_UNITS` (now merged into `UNIT_MAP` in `gui/backend/units.py`; the separate map is gone), for example `entrainment_ratio [-]`, `m_dot_primary [kg/s]`, `P_c_star [Pa]`, `dP_tap [Pa]`, `dP_permanent [Pa]`, `UA [W/K]`, `NTU [-]`, `sigma_norm [-]`, `T_static [K]`.
  - **Patterns:** the junction ports `port_<i>_<quantity>` (e.g. `port_0_P [Pa]`, `port_1_m_dot [kg/s]`), the tee's hidden steps `step_<common|straight|branch>_<quantity>` (e.g. `step_branch_dPt [Pa]`), and the arrays' per-row lists `rows_<key>`, which take `<key>`'s unit (e.g. `rows_T_hot [K]`, `rows_m_dot [kg/s]`).
  - **Newly given a unit:** `bypass_fraction`, `eta_r`, `eta_r_extrapolated`, `Cd_own`, `rows_P_backside` and `rows_unresolved`, which had none anywhere.
  - **Unchanged:** the meta columns (`id`, `label`, `type`, `kind`, the solve status, ...), the text-valued keys (`area_source`, `port_<i>_area_source`, `second_law`, `tap`, ...), and `sweep(..., metrics=[...])`, whose columns are the metric strings passed in. A test now exports one network per palette kind and fails on any column with no unit.

- **GUI node cards show curated rows per element kind (#775).** A card used to show the first Display Fields keys its result carried, in discovery order, so a channel card read `m_dot`, `Dh`, `M_exit_duct`. Each kind now has a short default list, most useful first: for example channel `m_dot`, `dP`, `mach_out`; orifice `m_dot`, `dP`, `Cd`; plenum `P`, `Tt`. The list lives in `gui/frontend/src/utils/cardRows.ts`, and a test checks every key against the backend's per-kind key list.
  - **Changing them:** a pin next to each Live Telemetry field adds it to (or removes it from) the cards of that kind. A new pin goes first, so it shows on a one-row card. "Reset card rows" restores the default.
  - **Fits the card:** the discrete loss card shows two rows (`m_dot`, `dPt`). A long Telemetry label gets a shorter card label (`ΔPt`, `ΔP`, `ṁs/ṁp`, `ṁp`, `ṁs`, `ṁbr/ṁcom`), and a test keeps every default row within the 21 characters a card shows.
  - **Not changed:** the Live Telemetry still shows every selected Display Field. The Display Fields panel no longer affects the cards.

- **The film density ratio takes the gas's static temperature (#782 audit).** `DR = rho_c / rho_g` is formed where the jets meet the gas. The gas moves at G_gas, so `rho_g` is at the static T of the duct closure at (T_gas, Pt, G_gas), as the effusion path's chamber density already was; it was at the gas Tt.
  - The coolant leaves a plenum and keeps its Tt.
  - Film effectiveness moves by under 0.7% at a gas Mach of 0.3 or less, and up to 2% at Mach 0.5. In the film test network, DR moves by -0.35%.
  - The scorecards are unchanged: none runs a film through `film_on_wall_side`.

- **Stagnation inputs are spelled `Tt`/`Pt` across the public API (#782).** The P/Pt rule of #773 (P static, Pt total) now holds for temperature: `T` is static, `Tt` total. Keywords `T0`/`P0` became `Tt`/`Pt` in `nozzle_flow`, `nozzle_cd`, `nozzle_thrust_cd`, `solve_A_eff_from_mdot`, `solve_P_back_from_mdot`, `critical_pressure_ratio`, `mach_from_pressure_ratio`, `mass_flux_isentropic`, `duct_static_state`, `duct_htc_static_state`, `duct_htc_static_state_dY`, `T_from_stagnation`, `P_from_stagnation` and the solver tools `sonic_static_pressure` and `orifice_compressible_mdot_and_jacobian`; `compressible.nozzle_flow(T, P, ...)` takes `(Tt, Pt, ...)`; `heat_entropy_flux` takes `Tt_in`/`Tt_out`; `approach_static_pressure` takes `(Pt, Tt, ...)`; the orifice and channel residual kernels take `Tt_up`. Functions `T0_from_static`, `T0_from_static_v`, `P0_from_static`, `solve_P0_from_mdot` and `T0_/P0_from_static_and_jacobian_M` are now `Tt_from_static`, `Tt_from_static_v`, `Pt_from_static`, `solve_Pt_from_mdot` and `Tt_/Pt_from_static_and_jacobian_M`. Fields: `MixtureState.Tt`/`.Pt` (were `T_total`/`P_total`), `MassStream.Tt`/`.Pt` (were `T`/`P_total`), `StreamJacobian.d_Tt`/`.d_Pt`, `NozzleSolution.Tt`/`.Pt`, `DuctStaticState.dT_dTt` and the other `d*_dTt`/`d*_dPt` partials, `SonicStatic.dP_dTt`, `OrificeResult.d_mdot_dTt_up`, `ChannelSolverResult.d_dP_dTt_up`. Positional calls are unchanged and every result is bit-identical. The old Python spellings still work until 0.9.0 (see Deprecated).
- **C++:** the same struct fields are renamed in `solver::Stream`, `solver::StreamJacobian`, `MixtureState`, `NozzleSolution`, `DuctStaticState`, `SonicStatic`, `OrificeResult` and `ChannelResult`, without an alias; the renamed free functions keep `[[deprecated]]` inline forwarders until 0.9.0.
- **Kept:** `Pt0`/`Tt0` of the Fanno and area-change functions (the `t` says total, the `0` is the inlet station beside `Tt_exit`), and the `T_from_stagnation`/`P_from_stagnation` names. Network result keys and `NetworkMixtureState` are not part of this change.

- **Network results report the STATIC temperature under `T` (#782), as they already did for `P`.** The same key name changed meaning.
  - `"<node>.T"` and every element's `T_in`/`T_out` are the static temperature. Until 0.7 they held the stagnation temperature, which is now `"<node>.Tt"` (on every node) and `Tt_in`/`Tt_out`.
  - A momentum chamber's `"<chamber>.T"` is the static T its closure solved, formerly `"<chamber>.T_static"`.
  - The static T is regime-aware. A plenum (`P = Pt`) and an incompressible node report `T = Tt`. A moving compressible node's `T` lies on the isentrope through (`Tt`, `Pt`) at its static `P`.
  - Reading such a key on a `solve()` result warns once (`FutureWarning`) where its value differs from the `Tt` one.
  - Solutions are unchanged: no residual reads `T` (#825). The residual vector and Jacobian are bit-identical on all 560 stress networks.
  - `extract_complete_states()` now pairs the static `P` with the static `T`.
- **`NetworkMixtureState.T` is a computed property, the static temperature.** The constructor takes `incompressible=` for a node in that regime.
- **GUI:** node cards show `Tt` by default; `T`, `T_in` and `T_out` are labelled static.

- **Breaking (0.8.0): the CSV export and `NetworkResult.to_dataframe()` carry the case's inputs (#784),** so an export reads back into the full results table with `load_results(path)` alone. The result rows and their `<key> [<unit>]` columns (#774) are unchanged; after them come `input:<key> [<unit>]` columns (each entity's node or edge data, nested fields dotted: `input:L [m]`, `input:composition.source`), `initial_guess:<key>` columns, and on a new row of `type` `network` (id `_network`) the network-wide settings and the solver options as `solver:<key>`. An edge with inputs and no result row gets a row of `type` `edge`. A reader that iterates over every row now meets the `network` row: filter on `type`. Format: `docs/GUI_TECHNICAL.md`, "CSV Export".

- **Breaking (0.8.0): `NetworkRunner.sweep()` returns the results table, one row per case (#784).** It used to return a flat frame: with `metrics`, the override columns, `success`, `final_norm` and one column per metric string; without, `to_dataframe()` stacked per case with a `_sweep_index` column. It now returns `(kind, entity, quantity)` columns, each case's full inputs included, with the row labels taken from the index of `params`. Update readers: `sweep_df["combustor.Tt"]` becomes `sweep_df["outputs"]["combustor", "Tt"]`, `sweep_df["success"]` becomes `sweep_df["outputs"]["_solve", "success"]`, and an override column `sweep_df["air_inlet.Pt"]` becomes `sweep_df["inputs"]["air_inlet", "Pt"]`. `metrics` now selects which outputs to keep. `load_results()` reads an old full-detail sweep frame into the new shape.

- **An orifice's bore against its pipe has a rule; the silent wall-hole fallback is gone (#788). Breaking in 0.8.0.**
  - **Before:** a bore at or above the pipe it sat in silently became a wall hole (beta = 0), per direction, with only a warning on the forward side. That is how the stress `random` family, several test networks and the grid example wrote a lumped Cd*A restriction.
  - **Now:** a bore equal to its pipe (to `ORIFICE_BORE_PIPE_RTOL` = 1e-6 relative, on either side) is **no orifice**: zero loss, the row is `Pt_in - Pt_out = 0` (row kind "p"), the flow passes unchanged. It is how an orifice is switched off without editing the network. A bore larger than its pipe is refused with a `ValueError`.
  - **Migration:** write a lumped Cd*A restriction as `OrificeElement(..., in_pipe=False)` (GUI: untick "In a pipe"; JSON: `"in_pipe": false`). It has no pipe diameter, so beta = 0 in both directions and there is no pressure recovery. It is the old fallback, now stated instead of inferred. Its forward rows match the old fallback bit for bit. Two cases differ: an orifice that was in-pipe in one direction and a wall hole in the other is now a wall hole both ways, and a metering correlation (Reader-Harris/Gallagher, Stolz) evaluates its Cd at the plenum convention D = 1 m instead of D = 10 d.
  - **Branch selection and the under-expanded junction exit stay open under #788.**

- **Messages are ASCII** (#791): `bernoulli_v2` now says `(v2^2 < 0)`, the quarter-wave resonator's cutoff error says `kL = ... > pi`, the GUI backend's boundary-swap and unconnected-node errors use `--` for the dash, and the C++ material database's `al6061` notes read `<300degC`. Display units in `units_data.h`, the `_core` docstrings and the GUI's unit table keep their UTF-8 symbols. Every C++ and Python source is now checked by the `ascii-only` pre-commit hook with no exclusions.

- **Impingement jet holes take a crossflow-dependent Cd by default (#800).** A strong crossflow along the gap bends the jets over and lowers their discharge coefficient. `ImpingementPlateElement` and `ImpingementArray` gain `cd_crossflow`, also exposed in the GUI:
  - `'rogers_hersh_1976'` (the new default): Rogers & Hersh (1976, NASA CR-2681) measured single-hole fits. It applies to any plate and stays capped at the hole's own Cd below Gc/Gj ~1.5.
  - `'florschuetz_metzger_su_1982'`: CR-3630 Table 4.1, the authors' own law for the plates they tested, and the one behind Table D.1's flows.
  - `'none'`: the previous constant Cd.

  C++ `rogers_hersh_1976_cd` and `florschuetz_metzger_su_1982_cd` return the Cd with analytic derivatives. The law is written in terms of the ideal jet velocity, so the hole's flow goes as dP^1.07 under a strong crossflow instead of sqrt(dP).
  - **Unchanged:** every Florschuetz (1981) array is bit-identical, since its Gc/Gj stays at or below 0.77, so the Eq. 8 flow model and the Nu chain (MAE 6.88%) do not move.
  - **CR-3630 Table D.1 (cross-source law):**
    - The stagnant upstream rows of B(5,4,1) at mc/mj 0.99 no longer alternate in sign (#799); they now rise monotonically, as the data do.
    - Worst-row Gj/Gj_mean improves: z/d >= 2 max 0.111 -> 0.085, z/d = 1 max 0.269 -> 0.231.
    - P/P0 mean error improves 0.0075 -> 0.0073.
  - **CR-3630 Table D.2:**
    - eta is unchanged (bias +0.013, RMSE 0.077).
    - Nu is unchanged inside the correlation's range.
    - Starved rows go from 8 to 11, with Nu at -81% (-95% before): there is still no channel-flow floor (#467).
  - **Solver:** a cold solve that fails on such a network is retried once from the constant-Cd root. CR-3630 D.2 B(5,4,1) at mc/mj 1.01 with a heated bypass needs it.

- **Crossflow-segment stations take the axial momentum of the jets that join them (#813).** A merge chamber already credited each side stream with m w cos(theta), using the injecting element's `injection_angle_deg` (an `EffusionPlateElement`'s `angle_deg`). A `CrossflowSegmentElement` station ignored the hole's angle and took the side stream's axial velocity as `kappa` times the receiving duct's own. Both now read the same declaration: an inclined hole joining a crossflow chain pushes the crossflow, ejector-like, instead of entering normal to it.
  - **Unchanged:** normal jets (90 deg, the default, and every impingement row) declare no momentum, so every existing `ImpingementArray` and `EffusionLiner` solves bit-identically, including the Florschuetz and CR-3630 scorecards.
  - **Second law:** a hole whose jet brings axial momentum is grouped with the chain's control volume, entering at its supply's total pressure. Closed at the station's Pt instead, the jet's momentum read as pressure the chain created, and a jet-pumped chain was demoted to inconsistent.

- **The combustor head loss is total-total, as in Lefebvre's pressure-loss factor (#782 audit).**
  - `ConstantHeadLoss`/`LinearThetaHeadLoss` evaluate the inlet dynamic head at the inlet stagnation state (`Tt_in`, `Pt_in`) and divide by `Pt_in`. They used the inlet's static P with Tt, which overstated the loss at a moving inlet: +8% at M 0.3, +24% at M 0.5.
  - At a plenum inlet (P = Pt) nothing changes: residual vector and Jacobian bit-identical on all 560 stress networks; scorecards byte-identical.
  - The loss context is a `PressureLossContext` carrying `Tt_in`/`Pt_in`. `combustor_head_loss` takes `Tt_in`/`Pt_in`, and `CombustorLossResult` reports `dxi_dTt_in`/`dxi_dPt_in`.

- **Copyrighted third-party documents live in a top-level `sources/`, never in git (#846).** Papers, handbooks, scans, figure crops and supplementary archives used to sit in gitignored folders under `docs/`, protected by a per-folder and per-extension denylist that leaked for every new file type or folder.
  - `sources/` is ignored whatever the file type; `sources/.gitignore` is a second lock.
  - `docs/` is an allowlist: only directories and `.md` files are un-ignored, so a new file type or folder there stays out of git until someone allows it.
  - A `no-sources` pre-commit hook and a CI step (`scripts/check_no_sources.py`) reject any committed `sources/` path and any document, archive, office or image file outside an allowlist of directories holding images we made. This catches `git add -f`, `--no-verify` and web edits.
  - `validation/sources.yaml` entries point at `sources/<topic>/...`; `scripts/sources.py check` takes `--sources` (was `--docs`) and reports register entries outside `sources/`, unregistered documents and any file there git would track.
  - `scripts/migrate_docs_to_sources.py` moves an existing checkout's untracked `docs/` files to `sources/` (dry run unless `--apply`; `docs/sources/` scans become `sources/scans/`). Run it once after pulling this change: an untracked `.md` left in a `docs/` subfolder is no longer ignored.

- **GUI: a probe keeps the solve, and shows its target by a dotted link (#856).** A probe is display-only (the backend drops it before solving), so adding, removing, pasting or retargeting one no longer clears every node and wall result and the solve status; a probe dropped onto a solved network reads its target as soon as it is pointed at one. Drag from the probe's right-hand dot onto a node or a thermal wall (its path or label box) to attach it; the attachment is drawn as a thin dotted, muted link to the node's border or the wall's label. Select the link and press Delete to detach the probe; picking a target in the Inspector moves it. The link is derived from the probe's `target_id`, so saved files are unchanged. Selecting or moving a node no longer clears the solve status either, which used to take the Inspector's second-law and Mach flags with it.

- **A sharp-edged orifice now expands like one: `jobson1955` is the default compressible model for incompressible-referenced Cd correlations (#861). Breaking in 0.8.0.**
  - **Before:** every `OrificeElement` used the isentropic nozzle unless `compressible_model` was set. A sharp plate's jet widens as the pressure ratio falls; the nozzle's flow is flat past choke and 3-9% low at 10-20% pressure drop.
  - **Now:** `compressible_model` defaults to `None`, and the element selects a model ONCE at build from its correlation and geometry, never per iterate. `"jobson1955"` for a Cd referenced to incompressible flow (`ReaderHarrisGallagher`, `IdelchikThick`/`Beveled`/`Rounded`, `Lichtarowicz`, `McGreehanSchotsch`) on a sharp-edged bore (`edge_radius = 0`, `bevel_depth = 0`) at beta <= `JOBSON_1955_BETA_MAX` (0.75) in both flow directions. `"nozzle"` otherwise: `"fixed"` (the measurement basis of your Cd is unknown), `"Stolz"` (ISO 5167:1980 not held), a rounded or beveled edge, a wider bore. The cooling-plate and connection subclasses (`EffusionPlateElement`, `FilmCoolingElement`, `SlotFilmElement`, `ImpingementPlateElement`, `EffectiveAreaConnectionElement`, `DiameterDischargeCoefficientConnectionElement`) keep the nozzle; McGreehan-Schotsch's crossflow term is Rohde-fitted and nozzle-referenced. A Ci above Jobson's 0.7 is flagged by `model_in_range`, not refused. The incompressible regime is unaffected.
  - **Which results move:** compressible networks with a plain `OrificeElement` on one of those correlations, including the default `ReaderHarrisGallagher`. Its flow rises by about the gap the nozzle left: through a default element, Perry's air (Jobson 1955 Fig. 5a) goes from RMS 20.8% to 1.3% and Schiller's steam (Fig. 5b) from 20.0% to 2.1% (cross-source). The #852 ISO 5167 metering example reads its own flow back within 0.1% instead of 0.8-0.9% high. Networks whose orifices keep the nozzle are bit-identical.
  - **To keep the old numbers:** pass `compressible_model="nozzle"`. A `Cv` other than 1 now needs `compressible_model="jobson1955"` set explicitly.
  - **GUI:** the Inspector's "Compressible Model" gains "Default (selected at build)" and shows what it selected. A saved file without the key, or with `"default"`, gets the new selection; a file that stored an explicit `"nozzle"`, `"jobson1955"` or `"iso5167"` keeps it.
  - **New constant:** `JOBSON_1955_BETA_MAX` (0.75), from C++'s `orifice::jobson1955::beta_max`.

- **Stall handover catches far-from-root crawls that cannot finish in time (#864).** A cold hybr or LM-fallback phase more than 100x tol from the root now also hands over to the next rung when the time to tolerance projected from its recent rate (over the last 5 s) is more than 3x the time left in the attempt, or 3 x 60 s without a timeout. Before, such a crawl improving just over 10 % per 5 s ran to the timeout. Example: `jet_chain` 32 under `jobson1955` failed at 20 s; it now converges in about 6 s through the warm-start retry.
  - **Measured:** serial runs on 14 stress families x 80 seeds and the 96-case audit set, under both orifice models. No success was lost, one was gained, and total time fell by 2-3 %.
  - **Opt-out:** `solve(stall_handover=False)` or the GUI checkbox turns it off with the rest of the handover.

- GUI: a node's position is now its card centre, and the canvas snaps to a 10 px grid. Cards of different heights placed on one row now line up their flow handles, so a main path runs straight. Saved networks are written as file version 1.2. Older files and the autosave are converted on load and land where they were drawn; the bundled examples are re-laid on one centre line (#865).

- **`Cd` on an orifice-family element is the `"fixed"` correlation's value only; a `Cd` given with any other correlation raises (#872). Breaking in 0.8.0.**
  - **Before:** `OrificeElement`, `EffusionPlateElement`, `FilmCoolingElement`, `SlotFilmElement`, `ImpingementPlateElement` and `ImpingementArray` took `Cd=0.6` (0.79 on the impingement elements) by default and, under any correlation but `"fixed"`, used it only as the solver's start guess. `OrificeElement(..., Cd=0.8)` ran the default `ReaderHarrisGallagher` Cd, about 0.6, with no warning.
  - **Now:** `Cd` defaults to `None`. `Cd` with a correlation other than `"fixed"` raises `ValueError` naming `correlation="fixed"`; `"fixed"` without a `Cd` raises too, except on `ImpingementPlateElement`/`ImpingementArray`, whose documented `"fixed"` default stays Florschuetz's 0.79. The connection elements (`EffectiveAreaConnectionElement`, `DiameterDischargeCoefficientConnectionElement`) are `"fixed"` with a Cd already and are unaffected.
  - **To migrate:** for a measured or assumed Cd add `correlation="fixed"`; to use the correlation drop `Cd`. Dropping it gives the numbers you had before. Adding `"fixed"` gives the Cd you wrote, so the results move.
  - **Start guess:** with a correlation it still assumes Cd 0.6 (0.79 on an impingement plate), so a network that never passed a `Cd` with a correlation solves bit-identically.
  - **Examples:** every example that passed a `Cd` with the default correlation meant a fixed Cd, and now says `correlation="fixed"`. Their printed numbers move:
    - `flow_regime_comparison.py`: the network's inlet pressure now matches the direct `cb.orifice_mdot` calculation it is compared with to 0.02-0.11%. It was 2.8-7.8% (incompressible) and 1.2-24.7% (compressible) off.
    - `combustor_network_example.py`: combustor Pt 130.9 -> 120.5 kPa.
    - `combustor_pressure_loss_example.py`: burner Pt 4.951 -> 4.846 bar.
    - `axial_staging_combustor_network_example.py`: the zones' Tt 2003.69 -> 2001.82 K.
    - `network_from_grid.py`: the largest grid's inlet pressure 1.5 -> 1.3 bar.
    - `combustor_liner_cooling_network.py`: its near-lossless exits (Cd 0.98 at 95% of the pipe bore) are also `in_pipe=False` now, lumped restrictions. In the pipe, a Cd of 0.98 at beta 0.95 left the coupled solve unconverged. Its printed results are unchanged; the coupled case's coolant dP rises 0.08%.
    - The benchmarks and the `docs/API_PYTHON.md` snippets are updated the same way.
  - **GUI:** a node keeps its `Cd` field whatever its correlation, and the Inspector shows it for "Manual / Fixed Value" only. The backend now passes `Cd` for `"fixed"` alone. A saved file with both a `Cd` and a correlation still loads and solves as before; its `Cd` is not used, as it never was.

- **Reader-Harris/Gallagher Cd is held below its lower Re_D bound instead of climbing to 1.5 (#873). Behaviour change in 0.8.0.**
  - **Before:** `metering_cd_and_dRe` evaluated ISO 5167-2 Eq. (4) far below its Re_D range. Its Re terms grow without bound as Re_D falls, so Cd rose until a 1.5 cap, with a Re_D < 1 floor under it. Both were derivative kinks. Example, D 4 mm, beta 0.5, D-D/2 taps: 1.5 at Re_D 1-40, 1.162 at 100, 0.687 at 1000. A near-zero-flow plain orifice ran at A_eff = 1.5 A.
  - **Now:** below the bound Re_b = `metering_cd_re_min(ReaderHarrisGallagher, geom)`, Cd is held at the correlation's value near Re_b. Over ln(Re_D/Re_b) in [-0.1, 0] the slope dCd/d ln Re_D ramps linearly to 0, so Cd and dCd/dRe_D stay continuous (C1) and Cd is monotone. The held value is C(Re_b) - s w/2, with w = 0.1 and s the log-slope at Re_b. It sits within half of ISO 5167-2 5.3.3.1's stated uncertainty of C, for every in-range geometry; that bound is what fixed w. The same example now reads 0.6345 everywhere below Re_D 4520. Another: D 100 mm, beta 0.6 goes from 0.728 at Re_D 1000 to 0.640.
  - **Bounds:** ISO 5167-2:2003 5.3.1: 5000, or 16000 beta^2 above beta 0.56; flange taps take max(5000, 170 beta^2 D[mm]). These are the ones `metering_cd_in_range` already used, now defined once; it still returns False below the bound, since the hold is an out-of-range placeholder (real sharp orifices lose Cd in laminar flow).
  - **Stolz is unchanged, bit for bit**, Re_D < 1 floor and 1.5 cap included. ISO 5167:1980's lower Re_D bound is not held as a source, so Stolz gets no hold and `metering_cd_re_min(Stolz, ...)` raises `ValueError` (C++: `std::invalid_argument`). Stolz is handled in #876.
  - **Unchanged:** every RHG Cd at Re_D >= Re_b, bit for bit. ISO 5167-2 Annex A (about 4800 entries) still matches to rounding. The F/J hashes of the 560 stress networks are identical with their own (fixed-Cd) orifices. With plain orifices forced to RHG, 130 differ, every one with a Cd call below the bound.
  - **Removed:** RHG's Re_D < 1 floor, which can no longer be reached. The 1.5 cap stays for geometry far outside the standard's: flange taps in pipes below D ~32 mm.
  - **New:** `metering_cd_re_min` (C++ and Python; Reader-Harris/Gallagher only), plus the constant `orifice::metering_hold::blend_ln_re`.
  - **Which results move:** a `ReaderHarrisGallagher` orifice running below its bound. A plenum-fed or `in_pipe=False` one evaluates Re_D on the D = 1 m stand-in pipe, so it is usually below. There, Cd is now close to the high-Re value and no longer follows the flow.

### Deprecated

- **`k_aluminum_6061` and the `al6061` material are deprecated, to be removed in 0.9.0** (#761).
  - Its 167 W/(m K) is the ASM Handbook Vol. 2 room-temperature value for T6 (O temper 180, T4 154, both at 25 C), but the falling temperature dependence has no source. The NIST cryogenic fit rises up to 300 K and gives 155 there.
  - Both names now emit a `DeprecationWarning` and still return the old values until removal.
  - Give the wall layer its conductivity directly (`material="custom"`) instead.

- **More stagnation quantities still spelled the pre-#782 way are renamed, with the old spellings deprecated until 0.9.0.** Each old keyword or field still works and returns the replacement's result bit for bit, with a `DeprecationWarning` that names its replacement; passing an old keyword together with its replacement raises `TypeError`.
  - `chamber_merge_face_state(T_main=, T=)` -> `Tt_main=`, `Tt=` (both always were stagnation temperatures; the docstring called `T_main` static).
  - Internal solver tools: `momentum_chamber_residual_and_jacobian(P_total=, T=)` -> `Pt=`, `Tt=`; `orifice_tap_pressure_and_jacobian(P_total_up=, P_total_down=)` -> `Pt_up=`, `Pt_down=`; the `p0=`/`t0=` keywords of `ejector_choked_mass_flow_and_jacobian`, `ejector_cd_nozzle_mass_flow_and_jacobian` and `ejector_cd_nozzle_exit_static_and_jacobian` -> `Pt=`, `Tt=`.
  - `MixerResult.T_mix`, `.P_total_mix`, `.dT_mix_d_stream`, `.dP_total_mix_d_stream`, `.dT_mix_d_delta_h`, `.dT_mix_dQ` -> `Tt_mix`, `Pt_mix`, `dTt_mix_d_stream`, `dPt_mix_d_stream`, `dTt_mix_d_delta_h`, `dTt_mix_dQ` (the mixed TOTAL state).
  - C++ only, renamed without an alias: `ChamberStaticState::lnP0P` -> `lnPtP` (and `ChamberStaticStateDX::dlnP0P_dX` -> `dlnPtP_dX`); the `MixerResult` fields above.

- **A residual Jacobian keyed `"<node>.T"` is deprecated; the node's stagnation-T column is `"<node>.Tt"` (#782).**
  - Every built-in element now uses `"<node>.Tt"`. The residual vector and Jacobian are bit-identical on all 560 stress networks.
  - A custom element that returns the 0.7 name `"<node>.T"` gets the same Jacobian, with one `DeprecationWarning` per process. The old name is removed in 0.9.0.

- **The `T0`/`P0`, `T_total`/`P_total` and stagnation `T`/`T_up` spellings renamed in #782 are deprecated and will be removed in 0.9.0.** Each old keyword, function name and struct field still works and returns the replacement's result bit for bit, with a `DeprecationWarning` that names its replacement; passing an old keyword together with its replacement raises `TypeError`. `docs/API_PYTHON.md` ("Stagnation names") lists them all; `python/tests/test_deprecations.py` checks each one.

- **`"<chamber>.T_static"`, `NetworkMixtureState(T=...)` (and the positional `T` of `(P, Pt, T, Tt, m_dot, Y)`), and assigning `state.T` are deprecated and will be removed in 0.9.0 (#782).**
  - `"<chamber>.T_static"` returns the same value as `"<chamber>.T"`.
  - The constructor argument and the assignment set `Tt`, which is what `T` held until 0.7.
  - Each warns with a `DeprecationWarning`.
  - The 0.8 `FutureWarning` on reading a changed `T` key goes in 0.9.0 too.

- **`PressureLossContext.T_in`/`.P_in`, `combustor_head_loss(T_in=, P_in=)` and `CombustorLossResult.dxi_dT_in`/`dxi_dP_in` are deprecated and will be removed in 0.9.0 (#782 audit).** Each reads or forwards to its `Tt_in`/`Pt_in` replacement with a `DeprecationWarning`. A hand-built 0.7 context with only `T_in`/`P_in` still evaluates the built-in correlations; note that its `P_in` is now taken as the inlet TOTAL pressure.

- **The Stolz metering discharge coefficient (ISO 5167:1980) is deprecated and will be removed in 0.9.0 (#876).** ISO replaced Stolz with Reader-Harris/Gallagher in ISO 5167-1:1991/Amd 1:1998 (and ISO 5167-2:2003, 2022). Reader-Harris/Gallagher has a sourced range of use and a bounded low-Re hold (#873), and fits the EEC/API database better over Stolz's own range (Reader-Harris & Sattary 1996, Table 10: s 0.245 % against 0.390 %). Stolz's limits of use are not held, so `metering_cd_in_range` cannot flag a plate outside them. Results are unchanged until removal.
  - `OrificeElement(correlation="Stolz")` emits a `DeprecationWarning` naming `'ReaderHarrisGallagher'`, once per element at construction (not per residual). The cooling plates that refuse the whole ISO 5167 family refuse Stolz without the warning.
  - `cb.MeteringCdCorrelation.Stolz` passed to `metering_cd_and_dRe`, `metering_cd_in_range`, `solve_orifice_mdot`, `orifice_flow` or `orifice_flow_state` emits a `DeprecationWarning`.
  - C++: `orifice::Cd_Stolz` is `[[deprecated]]`; `MeteringCdCorrelation::Stolz` is marked deprecated in comments. The orifice example dropped its Stolz column.
  - GUI: the option reads "Stolz (ISO 5167:1980, deprecated)" with a note under it. A saved network with a Stolz orifice still loads and solves, and the solver status shows a deprecation note for it (`SolveResponse.deprecations`).
  - Migration: use `correlation="ReaderHarrisGallagher"` (the `OrificeElement` default) or `MeteringCdCorrelation.ReaderHarrisGallagher`.

### Removed

- **`standard_dry_air_composition()`** is removed (Python and the `_core` binding). Use `species.dry_air()`.
- **`smr_wgs_equilibrium` / `smr_wgs_equilibrium_adiabatic`** are removed: Python, the bindings, the C++ functions and their `units_data.h` entries. Use `reforming_equilibrium` / `reforming_equilibrium_adiabatic`, which agree with them on CH4-only input; the C++ and Python tests now run those cases against the replacements, with the removed solver's values pinned as golden data.
- **`combaero.network.ConstantLoss` / `LinearThetaLoss`** are removed, with the module `__getattr__` that warned on them. Use `ConstantFractionLoss` / `LinearThetaFractionLoss`.
- **`BorderCarnotLossElement`** is removed, with its C++ backing (`border_carnot_loss_residual_and_jacobian`, `BorderCarnotLossResult`, and `include/multi_port_chamber.h` with `border_carnot_L` and `HAGER_FRACTION`). Use `PressureLossElement(id, from_node, to_node, correlation=ConstantHeadLoss(zeta=4 * (1 - cos(0.75 * delta)) ** 2, area=area), area=area)`, with `delta` in radians, which reproduced it to round-off. Never put one on a `MultiPortChamberElement` or `ConstantKTeeElement` port (#272). The GUI never offered it.
- **`VortexElement(omega_rpm=...)` and `VortexElement.omega_rpm`** are removed. Use `shaft_speed_rpm` (#773). A GUI file saved with `omega_rpm` still loads.
- **`VortexElement(r_c=..., n=...)`** are removed: the free vortex (#661) ignored them. Passing either is now a `TypeError`. The `cb.vatistas_*` functions are unchanged.
- **The Mynard / constant-K tee's old `initial_guess` key `{id}.P_jct`** no longer seeds the junction total. Use `{id}.Pt_jct` (#773). `legacy_unknown_names()` is gone. A saved GUI file's tee guess is still moved to `Pt_jct` on load.
- **`init_strategy='incompressible_warmstart'`** is removed and raises `ValueError` listing the valid strategies. Use `'outletref_warmstart'`, or leave the default. Its only knob, `NetworkSolver.solve(warmstart_maxfev=...)`, went with it. The GUI backend maps a saved `incompressible_warmstart` to `outletref_warmstart`, as the frontend already did.

### Fixed

- **Rich equilibrium no longer stalls once the steam is used up (#759).** `combustion_equilibrium`, `reforming_equilibrium(_adiabatic)` and `combustion_equilibrium_and_jacobian` solve the steam-reforming extent of each hydrocarbon and the water-gas shift jointly -- Newton on ln Q - ln K = 0 with the analytic Jacobian, rounds of single-reaction solves only as the start, an element-potential restart (Gordon and McBride, NASA RP-1311) where Newton makes no headway -- instead of repeating rounds until nothing moves (#755).
  - **Before:** rich propane at phi >= 2.5 burned to a mixture with no steam left; the rounds got it back through reverse shift a sliver per round and stopped after 2000, `converged=False`, 124 / 229 / 320 K above Cantera's HP equilibrium at phi 2.5 / 2.75 / 3.0 (300 K inlet), and the derivatives chained the unconverged rounds (dT/dT_in 1.5 % off).
  - **Now:** converged, within 1.5 K of Cantera (gri30) and 1e-4 in mole fraction for CH4 and C3H8 over phi 1.5-3.0; the derivatives are one implicit-function solve on the converged system and match central differences to 3e-7. Unchanged elsewhere: the #755 reforming cases move by under 3e-13 in mole fraction and 1.2e-10 K, CH4 and lean/moderate propane by round-off. The steam-starved solve takes 1-4 ms instead of 0.2-0.3 s.

- **`k_haynes230` follows the manufacturer's table (#761).** Its conductivity now matches the HAYNES 230 brochure H-3000O: 8.9 W/m K at room temperature, rising 2.0 per 100 C to 28.4 at 1000 C, which is exactly k = 8.4 + 0.020 T[C]. The quadratic it replaced was 37% high at room temperature and 2-6% high at liner metal temperatures (900-1200 K). The valid range now ends where the table does, at 1000 C (1273 K) instead of 1400 K: beyond it `k_haynes230` raises and the database lookup holds 28.4.

- **GUI: a cleared Max Solver Time is no limit, and a timed-out solve no longer seeds continuation (#762).**
  - The empty Max Solver Time field now reads "No limit" (it read "None").
  - `timeout: null` on `/solve` and `/export` now sets no request deadline either; the route used to substitute 180 s and answer 504 at 190 s. `timeout: 0` bounds the request at what the solver itself gets (0 s, it times out at once) plus the 10 s grace, instead of 180 + 10 s.
  - The 504 now names the Max Solver Time and the grace it waited, and suggests raising it or clearing it for no limit. It is only sent when a limit was set.
  - A worker thread that finishes after its request answered 504 can no longer write the continuation state: each request holds a token while it waits, and the state is stored only while that token is live. The thread itself still cannot be cancelled and runs to completion in the background.

- **NaN input is refused instead of answered (#769).** Domain checks written `if (x <= 0) throw` let NaN through, because every comparison with NaN is false: the function returned NaN, or after a clamp a plausible number. They now use NaN-safe predicates from a new public header, `input_validation.h` (`combaero::validation::is_finite_positive`, `require_finite_positive<E>`, ...), so NaN and +/-inf raise. Each function keeps the exception type it already raised (`ValueError`, `RuntimeError`, or `IndexError` for the saturation-pressure range).
  - **Humid air:** `saturation_vapor_pressure`, `vapor_pressure`, `humidity_ratio`, `water_vapor_mole_fraction`, `dewpoint`, `relative_humidity_from_dewpoint`, the functions built on them, and `HumidAir.set_TP_RH` / `set_TP_omega`. `humidity_ratio` and `water_vapor_mole_fraction` now check P (finite and > 0, `RuntimeError`).
  - **Materials:** the `k_*` functions; `get_material_conductivity(name, nan)` returned 1e-3 W/(m K) (clamped and floored) and now raises `RuntimeError`. A finite T outside the range is still clamped.
  - **Acoustics:** Helmholtz, tube and annulus modes, Q and layer thicknesses, `wavelength`, `strouhal`, `acoustic_properties`, `orifice_impedance_with_flow` (all arguments, including the bias and grazing velocities) and `is_whistling_risk`.
  - **Geometry:** `annular_area`, `hydraulic_diameter*`, `circular_area`, `cylinder_volume`, residence time and space velocity. The `Tube`, `Annulus` and `CanAnnularFlowGeometry` methods did not validate at all -- a swapped annulus returned a negative area -- and now raise `ValueError` unless D_outer > D_inner >= 0 and the sizes they read are finite and positive.
  - **The same pattern elsewhere:** existing checks in the incompressible, orifice, stagnation, compressible, friction, heat-transfer, cooling, combustion, thermo, transport, state and Vatistas functions. +inf stays accepted where it is a physical limit (NTU, a film coefficient h, YSZ aging hours). Converged network solves are unchanged.

- **The compressible orifice's back-pressure derivative near choke (#772).** `orifice_compressible_mdot_and_jacobian` scaled d(mdot)/dP_back by a cubic blend over 0.01 in pressure ratio either side of the critical ratio, but returned the unblended flow, so in that band the Jacobian was not the derivative of the value.
  - **Before:** d(mdot)/dP_back up to 46% off a central difference of mdot at P* + 50 Pa (the blend weight reached 1.5 at the band's upper edge and dropped back to 1 there), and d(mdot)/dPt, which reuses it, up to 0.3% off.
  - **Now:** the exact derivative of the value as written: (M^2 - 1)/v on the isentrope when unchoked, zero when choked. The flow is already C1 through P* (the mass flux peaks there), so no blend is needed. Flows are unchanged; converged networks agree to round-off.

- **Two reported Reynolds numbers now take the viscosity at the local static temperature.** Diagnostics only: no flow, loss or solver result changes.
  - `fanno_channel` (constant f) stored each profile station's `Re` with the inlet viscosity: 1.9 % off at the outlet of an air duct from 900 K, M 0.5 to M 0.66. The rough overload already used the local one.
  - `area_change_loss`, conical branch: `Re` used mu at the stagnation `Tt0`, 4 % low at M 0.62; it now uses the closure station's static T, as the sharp branch and `area_change_flow` do.

- **Kernels documented as taking a static T now get one (#782).** A network node carries its stagnation temperature, and four kernels that form a density from the static P beside them were handed it: the compressible tee's `BranchInput.T` (`TeeJunctionElement`), `station_half_drop` and `channel_entry_drop` (`CrossflowSegmentElement`, used by `EffusionLiner` and `ImpingementArray`), and `crossflow_velocity_ratio` (a duct-fed `EffusionPlateElement`).
  - **Before:** each density was low by Tt/T - 1: about 2 % at M = 0.3, 5-10 % at M ~ 0.5-0.6.
  - **Now:** each gets the static T of its node's moving gas, from the momentum chamber closure at the arm's or the station's mean flow (a plenum or boundary port is at rest: static = total), consistent with the kernels' documented static inputs. The Jacobians chain the static T onto the node's T0, static P, composition and flows analytically.
  - **CR-3630 Table D.1 (fidelity):** the worst-row Gj/Gj_mean error at z/d >= 2 moves from median 0.0205 to 0.0212 and max 0.100 to 0.111, inside the report's +/-2 % measurement band on Gj/Gj_mean (0.04 at that row); P/P0 mean error unchanged at 0.007.
  - **B(5,4,1) mc/mj=0.99 (D.1) is no longer reported as solved:** there was no root before -- the reported |F| 1.9e-4 was the overshoot past a false Fanno choke at the exit segment (#802: a station's static P is fed to the kernel as P0). With this change a regular root exists, but a cold start still stalls on that choke, as it now does for the D.2 test B(5,4,1) mc/mj=1.01. Both are strict xfails until #802 (#799).

- **Crossflow stations and area combustors report their moving state (#782 audit).** These nodes hold Pt = P while their gas moves, so they reported `.T` = Tt, `.Pt` = static P and `.mach` = 0. In one impingement network, a station at Mach 0.42 reported a Pt 12.8% low and a T 20 K high.
  - A station now reports the static T and Mach of the station's mean flow, the state its half-drops use, and its true Pt.
  - The crossflow segment leaving or entering a station reports that Pt, together with `dPt` and `pr_total`. A plate's or hole's jets discharging into the station keep their drop to its static pressure, the pressure the jets see.
  - A combustor with an `area` reports the static T at G = m / area.
  - Reporting only: residual vector and Jacobian bit-identical on all 560 stress networks; scorecards byte-identical.

- **Duct-fed effusion panels below McGreehan-Schotsch's Re floor now have an exact Jacobian (#760, #801).** An `EffusionPlateElement` (or any discharge-hole orifice) solves its own flow consistently, m* = Cd(Re(m*)) K, and its row's sensitivities divide by 1 - K dCd/dm. Below the correlation's Re floor (1e4 for McGreehan-Schotsch, 25 for Idelchik) the Cd value is held, but that denominator used the slope C++ continues past the floor as a Newton aid.
  - **Before:** every panel row of a duct-fed `EffusionLiner` was uniformly 0.987x the finite differences (1.26e-2 in the row-scaled metric, the GUI palette's effusion liner) and 0.97x at duct Mach 0.5 (2.98e-2), at every FD step.
  - **Now:** the sensitivities take the held value's zero slope, and the panel rows agree with finite differences to under 1e-8. Converged solutions are unchanged; only Newton's path can differ.
  - New `discharge_cd_re_floor(correlation)`: the hole Reynolds number below which the Cd value is held.

- **Crossflow segments run their Fanno duct from the true stagnation state, not a station's static pressure (#802).** `CrossflowSegmentElement`, used by `ImpingementArray` and `EffusionLiner`, fed a station node's pressure to the compressible Fanno kernel as the duct's P0. That pressure is the channel's static pressure, so P0 was a dynamic head too low and exit segments choked far too early.
  - **Before:** CR-3630 D.2 B(5,4,1) at mc/mj 1.01 ran its exit segment at 99.7% of that false choke (81% of the true one). D.1's B(5,4,1) at mc/mj 0.99 had no root, and both cases stalled from a cold start. A 6-row array that passes 0.041 kg/s raised out of `solve()` once its free exit pressure ran negative on the choke plateau. The diagnostics reported every exit segment choked (`M_exit_duct` 1.0) even at M ~ 0.35.
  - **Now:** the station half-drops move the duct's end STATIC pressures. The kernel runs between those statics' stagnation pressures at the segment's flux, from the momentum chamber closure, with analytic partials. A terminal end feeding the segment keeps its own Pt. A plenum or boundary the segment discharges into keeps static continuity, and an end below the sonic static is held on the subsonic branch. `choked`, `M_in_duct` and `M_exit_duct` report the duct the row solves.
  - **CR-3630 Table D.1 (fidelity):** no row of the 41 tests that solved before moves by more than 3e-4 in Gj/Gj_mean (band +/-2%) or 1e-4 in P/P0 (print floor 0.005). B(5,4,1) mc/mj 0.99 now solves cold and returns to the score. Its worst row is 0.269 (row 4, a stagnant row, #799), and the score over all 42 tests has P/P0 mean error 0.0075.
  - **CR-3630 Table D.2 (cross-source):** B(5,4,1) mc/mj 1.01 solves cold and returns to the score. eta over 300 rows has bias +0.013 and RMSE 0.078; the 290 rows without that test are unchanged at +0.010 / 0.075. Nu counts are back to 196 inside, 104 outside and 8 starved.
  - **Florschuetz 1981 Fig. 6 chain:** MAE 6.88% -> 6.88%. **Effusion liner** (no data, model only): inlet flow changes by -0.1% at duct M ~ 0.09, -1.0% at M ~ 0.3 and +0.8 to +1.3% at M ~ 0.6, where the bleed fraction falls from 7.2% to 6.2% (4 stations).
  - **Reversed out of an exit plenum or boundary, the incompressible regime now charges the same acceleration (#802).** A crossflow segment whose flow reverses out of an exit node at rest drew it as a reservoir in the compressible regime but paid friction only in the incompressible one, so the incompressible reversed flow ran about 5x the compressible one. It now pays one dynamic head (no entry loss), and the regimes agree within 0.5% on a 1 kPa reversed draw.

- `NetworkSolver.solve()` no longer raises out of a solve that ends on a node at or below zero pressure (#810). A root the residuals close with a node pressure at or below zero (seen on a crossflow array's Fanno continuation, the free exit at -113 kPa) is now reported INCONSISTENT and listed in the new `__P_nonpositive__`, before any diagnostics evaluate it. Previously, element diagnostics raised `ValueError: pressure must be positive` there before the second-law check could demote the root. On a failed solve, a diagnostics error is recorded in `__caught_exceptions__` (phase `diagnostics`) instead of being dropped silently; on a successful solve it still raises.

- **Crossflow stations hold for reversed flow (#815).** `station_half_drop` used m|m| for the channel's momentum flux. The flux, m u = m^2/(rho A), is positive whichever way the gas moves, so the drop now uses m^2. Forward flow is unchanged.
  - **Before, with reversed flow** through a station:
    - a merge in a reversed chain RAISED the static pressure along the flow, acting as a pump: a mirrored chain ran 1.7x the flow, with its stations below the exit pressure;
    - a station the flow leaves both ways (a side stream splitting) had no physical root. Incompressible did not converge, and compressible landed on a root the second law rejected.
  - **Now:**
    - a reversed chain is its forward mirror: identical compressible, within 0.5% incompressible;
    - a split station solves symmetrically and passes the second law.
  - **Validation:** every scorecard is unchanged (CR-3630 D.1 and D.2, Florschuetz 1981), since their channel flows run forward.

- **`jet_impulse` returns analytic partials (#816).** The jet's streamwise momentum J = m w, read by merge chambers (#471), crossflow stations (#814) and duct-fed holes' U1/Vi (`crossflow_velocity_ratio`), took its dJ/dPt, dJ/dTt and dJ/dP from central differences inside C++.
  - **Before:** choked, they were about 1% off. The sonic point came from a golden-section mass-flux maximum, accurate only to about 1e-8 of Pt on its flat peak, and that noise sat inside a Pt*1e-7 step.
  - **Now:** they are exact on the variable-cp ideal-gas isentrope, unchoked and choked, and C1 across choking. The sonic point is located on that isentrope at M = 1 (T* from h(Tt) - h(T*) = a(T*)^2/2, P* = Pt exp((s0(T*) - s0(Tt))/R)); it agrees with `critical_pressure_ratio` to about 1e-8.
  - The values do not move beyond that last digit. A choked angled jet into a crossflow chain now passes a 1e-6 FD Jacobian check.

- **Element diagnostics are evaluated at the static state of their end nodes (#824).** An element's reported inlet/outlet `rho`, speed of sound, Mach and transport properties were formed at (`Tt`, static `P`). At a moving compressible node, T was up to 11% below Tt, so they were off by about that much.
  - A junction's `port_<i>_T` is now the port's static T; the stagnation T is the new `port_<i>_Tt`.
  - The ejector's `mach_jet` takes the primary's static speed of sound.
  - Solutions are unchanged: residual vector and Jacobian bit-identical on all 560 stress networks.

- **`orifice_flow` evaluates the density of the composition it is given (#829).** It formed rho from air's molar mass for every `X`, so its mass flow was off by sqrt(MW/MW_air) for any other gas: -26% for methane, for example. An empty `X` is dry air, as documented, now from the library's dry-air composition and gas constant instead of literals. That moves the air default by 0.02%.

- **A discharge-hole orifice now always solves its Cd(Re) flow to a root (#861).** `OrificeElement` solves m* = Cd(Re(m*)) |K| itself for a discharge-hole correlation (#582). Where Cd rises steeply with Re (Lichtarowicz, a long hole at low Re), the loop takes the slow fixed-point step, and its 50 iterations could end short of the root. That stop was returned silently.
  - **Before:** on the stress networks with Lichtarowicz forced on their orifices, 73 such stops were more than 1e-6 off, up to 1.7%. In one captured case the row's m* was 2.34e-5 kg/s against a root of 1.93e-5. No converged outcome changed in that sample.
  - **Now:** a stop at the cap that is no root is re-solved on a bracket, the same way as the in-pipe own-flow solve (#863). That solve can now also complete a one-sided bracket by probing, which the merged version could not.
  - **Unchanged:** every solve that met its root. The F and J hashes of all 560 stress networks are identical, the stress outcomes are identical, and no test result moved.

- **An orifice in a pipe now always solves its own flow to a root (#861).** `OrificeElement` solves m* = m_calc(m*, x) itself when its drive reads its own flow through the approach pipe (#651). Where m_calc(m) has a kink, Newton went back and forth across it until the 60-iteration cap, and the value it stopped at was not a root. The solver was never told.
  - **Where it showed:** two cases.
    - The approach pipe's choke at beta 0.95 (any `compressible_model`).
    - Jobson's held discriminant past its stated beta (`compressible_model="jobson1955"`, beta 0.84).
  - **Now:** when the cap is reached with |a - m_calc(a)| above 1e-9 relative, the root is solved by Brent's method on the bracket that Newton crossed. A solve that already met its root is unchanged bit for bit:
    - On the stress families the bracket solve never runs for the default nozzle.
    - The F and J hashes of all 560 stress networks are identical.
  - **Changed result: the GUI combustor -> channel -> orifice regression net** (`test_supersonic_chamber_root.py`; 33 kg/s through a 95 mm orifice in a 100 mm channel).
    - **Before:** the first attempt "converged" with the channel-to-orifice junction at Mach 1.46 and was refused as supersonic. The retry reported Mach 0.63 at junction Pt 56.7 bar. Both points carried an orifice row that was not a root: the orifice's consistent flow at 56.7 bar is about 38 kg/s, not 33.
    - **Now:** the solve converges directly to the network's consistent root, junction Pt 49.6 bar, Mach 0.86, and the combustor upstream at 52.1 bar instead of 58.3. The supersonic-root refusal itself is unchanged and is still exercised by the junction runner's Wang cases.
  - **Jobson stress bleed station (seeds 6 and 39):**
    - **Before:** the orifice row was 6-12% off and its Pt column 6% off a central difference. LM took 685 evaluations where the nozzle took 437, and both seeds ran out of their 20 s.
    - **Now:** both converge in 8.7 s.

- **A network solve no longer depends on the process's string-hash seed.** Three near-zero-flow blends (`OrificeElement` reverse flow, `PortDischargeElement`, the incompressible `ChannelElement` reversal) built their Jacobian dicts from a set union. A set of names iterates in hash order, which Python randomises per process, and the key order set how coinciding Jacobian terms were summed: a 1-ulp change that moved the Newton path. Affected networks reached the same root but differed from one run to the next within solver tolerance (up to about 4e-8 relative; final norms 1e-12 to 7e-4). The geometry-inference fallback that picked a neighbour's area from a set union is ordered too.

### Documentation

- **`k_tbc_ysz`'s docstring now describes the shipped model** (#761). It still cited the untraceable "NASA TM-2010-216765" and the replaced APS/EB-PVD behaviour. It now gives Zhu & Miller 2000 Eq. 1 and records the model's fidelity on the paper's own Fig. 4 at 32 h: it matches the 1320 C curve, and is 2% high at 1100 C and 3% high at 990 C. Eq. 1 is kept as printed.

## [0.7.1] - 2026-10-27

### Fixed

- **The Linux wheel builds again (manylinux_2_28, GCC 14).** GCC 14 at `-O3 -Werror` raised a false-positive `-Wstringop-overflow` on a vector allocation in `combustion_equilibrium_and_jacobian`. The Linux wheel failed, so **0.7.0 was tagged but never published to PyPI**, neither `combaero` nor `combaero-gui`; 0.7.1 is the first 0.7 release on PyPI and contains everything listed under 0.7.0. The Linux C++ CI job now builds with GCC 14, the wheel's compiler, so this is caught on a pull request.

## [0.7.0] - 2026-10-27

### Highlights

A long release (339 commits since 0.6.0); the sections below are the record.
What changes numbers or code on upgrade:

- **Compressible is the default regime (#585).** Networks built without a
  regime now solve compressible; pass `regime="incompressible"` per element
  (or set the GUI's global regime) to keep the 0.6 behaviour.
- **Transport and thermo properties move.** Air viscosity was 11.8% low
  (flagged NUIGMech Lennard-Jones parameters, #485); collision integrals are
  now exact and viscosity and conductivity C2 in temperature (#611). Species
  thermo defaults to the C2 re-fit of NASA-9 (#549); `"nasa9"` is still
  selectable.
- **Area changes use Idelchik's closures in a compressible, mass-flow form
  (#731).** The Mach correction and its `Mach` argument are removed.
- **Cooling correlations were rebuilt on provenanced parameter sets (#339).**
  Every correlation whose source did not survive review was removed, as were
  the unsourced orifice Cd arms (`'Miller'`, `'Auto'`, `'ThickPlate'`).
- **Deprecated names warn and are removed in 0.8.0**, among them
  `BorderCarnotLossElement` and `VortexElement(r_c=, n=)`.

New: film, effusion, slot-film and impingement cooling elements with scored
validation; ideal momentum tee; Sanger's incompressible jet pump and Huang
1999 ejector anchors; a Fanno-Rayleigh duct; a per-element second-law check
and network balance ledger; Newton-homotopy and boundary-condition
continuation retries; geometry inheritance and autosave in the GUI.

### Added

- **`OrificeElement(tap=...)` and the orifice's tap diagnostics (#645).** `tap` (`"Corner"`, `"D_D2"` or `"Flange"`, the default) picks the tapping differential the element reports as `dp_tap`; it never moves the solved flow. Between pipes the element also reports `dp_permanent`, `loss_ratio` and `recovery_in_range`, and `pressure_recovery` says whether the permanent-loss row applies.
- **`permanent_loss_ratio_and_dC(beta, C)`**: ISO 5167-2's permanent-loss ratio and its C derivative, and the solver's `orifice_tap_pressure_and_jacobian` (C++ and `combaero._solver_tools`).

- **An opt-in Newton-homotopy last retry for a cold solve that did not converge (#642):** `NetworkSolver.homotopy_retry = True`.
  - **How it works:** H(x, λ) = F(x) − (1 − λ) F(x0) from the cold guess (a state with flow), λ walked adaptively to a plain λ = 1 solve. Budget `NEWTON_HOMOTOPY_BUDGET` = 5 s.
  - **Who uses it:** the GUI backend's interactive solve and the stress harness turn it on. In the GUI it runs within the user's max solve time, all of what the cold solve leaves (`homotopy_budget = None`), so an unsolvable network reports once that time is up. The Sidebar's help text says so.
  - **Measured:** 3 of the 23 stress random draws that failed every other start now solve within 20 s.
  - **Why it is off by default:** a non-converging solve pays up to its budget.

- `PortDischargeElement`: a momentum chamber discharging into a `PressureBoundary`. Subsonic, its static pressure is the boundary's; choked, it sits at the sonic static pressure of its flux. Inflow, and a `coupling="total"` boundary, keep the total pressures equal.
- `_solver_tools.sonic_static_pressure(T0, G, Y)` (C++ `sonic_static_pressure`, `combaero::sonic_static_state`): the static pressure at which the mass flux G is sonic for stagnation T0, with variable cp and its derivatives.

- **The GUI's thermal wall declares its orientation (#629).** The wall
  Inspector gains an Orientation field (aligned or opposed) and shows each
  side's declared run, from the element's upstream to its downstream along
  the flow edges, never from which side a handle is drawn on. After a solve
  it shows the arrangement the flow signs gave (parallel or counter), with
  NTU and effectiveness. The wall edge carries an orientation marker, plus the
  solved arrangement. A saved file without the field stays aligned. Walls on
  an impingement array or effusion liner show no field: their row walls
  have no orientation to declare.

- **GUI autosave (#596).** The network (nodes, edges and solver settings; no solve results) is saved to the browser's local storage one second after the last edit, and at once when the tab is hidden or closed.
  - **Reloading the page** restores it, with a note saying when it was saved.
  - **The new Clear button** (two clicks) empties the canvas and deletes the saved copy.
  - **Storage is unavailable** (blocked site data, a private window, a full quota): the editor works as before, without the autosave.
  - **Verified in** Chromium, Google Chrome, Firefox and WebKit.

- **Shaped-hole film effectiveness (#386)**: `film_effectiveness_colban_2011` and its `(eta, deta/dM)` form, Colban, Thole & Bogard (2011) Eq. (19) for one row of fan or laidback fan holes. Baldauf (2002) stays the cylindrical-hole correlation; the two are not interchangeable at high blowing.
  - Lands on the Eq. (19) line the paper draws in Figs. 6-7 within 0.003 in eta (about a pixel). On those figures' Gritsch data, which the correlation was fitted on, it scores 3.6% RMS at M = 0.5 and 12.1% at M = 2.5.
  - Table 4's limits, plus the fitted angle (30 deg) and density ratio (1.7-2.0), are reported as extrapolation, not enforced.
  - The figures are committed as `validation/cooling/data/colban2011`, scored by the new `shaped_film_runner`.

- GUI: the Effusion Plate carries the same "f" film handle as Film Cooling Holes, at the same place (6 o'clock): drag it onto a thermal wall's label box to cover that wall with the plate's film (#386).

- **An effusion panel's film now covers the wall downstream of it (#386).** `ThermalWall(film="<panel id>")` accepts an `EffusionPlateElement` with `panel_length`.
  - Its rows are combined as on the panel's own wall: Baldauf (2002) for each row and Sellers superposition over them. The result is then read past the last row, averaged over the wall's span.
  - The panel's own wall is unchanged.

- **`FilmCoolingElement` and `ThermalWall(film=...)`: a film-cooled wall downstream of the holes (#386).**
  - The element is the holes: an orifice of `n_rows * holes_per_row` holes from the coolant supply into the gas, using Colban (2011) for shaped holes or Baldauf (2002) for cylindrical ones. It owns no wall.
  - A wall on the gas duct downstream names it with `film=`. That side's adiabatic wall temperature becomes `T_aw - eta_bar (T_aw - T_coolant)`, averaged over the wall's span (`film_start`, `film_end` from the last row). Its h, stream and heat booking are unchanged.
  - The Jacobian reaches the film's flow and the coolant temperature analytically.
  - The coupled wall is not scored; each kernel is checked on its own author's data.
- **`film_span_average_baldauf_2002` / `film_span_average_colban_2011`:** rows combined by Sellers superposition and averaged over a streamwise span, with exact derivatives in M and density ratio.
- **Every wall now reports `{wall}.T_aw_a` and `{wall}.T_aw_b`.**

- **GUI: Film Cooling Holes node and a thermal wall's Film (#386).**
  - The node is the holes, shaped (Colban) or cylindrical (Baldauf). It defaults to Colban's Gritsch configuration B, the geometry the correlation is scored on.
  - A thermal wall edge picks a film provider (film holes, or an effusion plate), the side it covers and the wall's span after the last row. Its results show η, T_aw before and after the film, M, the density ratio, and whether the film is outside its correlation's fitted range.

- GUI: a thermal wall's film shows as a dotted link from the film-cooling holes to the wall's label box while either is selected, and a "film" tag on the label box (its tooltip names the holes and the duct covered). Dragging from the holes' new "f" handle onto a wall's label box sets that wall's film. The film-side dropdown names the covered duct instead of "side A/B", and picking a film defaults it to the duct the holes feed (#386).

- An impingement array can take an initial (bypass) crossflow into its gap's upstream end: `ImpingementArray(..., bypass_length=None, bypass_entry_K=0.0).add_to(net, supply, exit, bypass=<node>)`, built on `ImpingementCrossflowElement(bypass_inlet=True)`. Fed from another array's exit, the gap continues with its momentum (arrays in series); fed from anything else, the bypass is accelerated from that node's total pressure. Row 1's Gc/Gj is the bypass over its jets. The flow split is checked against the authors' own data, Florschuetz, Metzger and Su (1982, NASA CR-3630) Table D.1, 42 tests: worst row per test median 0.021 in Gj/Gj_mean at z/d 2-3, up to 0.27 at z/d 1, where the authors' own constant-Cd model misses by up to 1.2. Heat transfer stays referenced to the supply, valid for a bypass at the jet temperature; the summary flags `heat_transfer_outside_source` otherwise. In the GUI the array node has a "b" handle; unwired, the end is closed (#467).

- Impingement arrays with a bypass at another temperature: each row's target reference temperature is `T_aw = T_j + eta_r phi_0 (T_0 - T_j)` (Florschuetz and Su 1985, NASA CR-3936, Eq. 2.11), with `phi_0` the bypass share of the crossflow approaching the row. `cb.impingement_eta_r(Gc_Gj, xn_d, yn_d, z_d)` is a closure fitted to CR-3936's inline Appendix A data, which publishes no correlation: `eta_r = X^n/(1+X^n)`, `X = a (Gc/Gj)(z/d)^r`, RMSE 0.083 leaving one geometry out. It is analytic, and the solver relays the reference temperature's dependence on the bypass flow and temperature into the Jacobian. Without a bypass, or with one at the supply temperature, results are unchanged. Diagnostics add `eta_r`, `bypass_fraction` and `eta_r_extrapolated`, and the array summary adds the rows' `T_aw`, `eta_r` and `bypass_fraction` (#467).

- **Ideal momentum junction kernel (#517)**: `_core.ideal_junction_residuals_and_jacobian`, a parameter-free junction closure with its whole-element (f, J). It is a baseline, parallel to the Mynard closure, and is not yet wired into a network element.
  - **Physics:** suppliers share one junction static pressure; each receiver takes its momentum along its own axis from the supplied mean velocity. A merge is Blaisdell & Manson's theory. A split loses the vector Borda-Carnot head on every outlet, which is Idelchik's diverging-branch Eq 7-6 with A' = 1. Nothing is fitted.
  - **Solver-friendly by construction:** one signed formulation covers merge, split and reversal, C-infinity everywhere including zero flow. There is no flow ratio, no direction classification and no degenerate-state fallback.
  - **Compressible:** exact momentum in static pressure with each port's own density. There is no Mach inversion in the kernel.

- `IdealMomentumTeeElement`: the parameter-free ideal momentum junction (#517) as a network element beside `MultiPortChamberElement` and `ConstantKTeeElement`. Its `P_jct` is the junction static pressure, seeded one common-port head below the common Pt. It takes 2 to 6 ports and has no declared flow direction. On 400 random-BC junctions within Mach 0.3 it converges on 200 of 203 (the shipped tee: 193 of 194), with a median solve about a third as long. A converged state with a reversed port is reported inconsistent, because junction port nodes cannot yet carry reversed flow.
- `validation.junction.random_robustness` takes a `model` argument (`mynard` or `ideal`); both models see the same draws.

- GUI: the tee node's junction model gains "Ideal momentum (parameter-free)" (`junction_model: "ideal_momentum"`), the #517 baseline. It needs no inputs beyond angle and areas. Flow direction only wires its ports, and the Inspector says so. Results show the junction static pressure `p_jn` and `n_reversed_ports`.
- `IdealMomentumTeeElement` diagnostics name a tee's flows by its geometry: `m_dot_com`, `m_dot_straight` and `m_dot_branch`, signed in the wired direction. K aliases appear only when the converged regime is the wired one, so a tee wired as a split that ran as a merge, because a fast common jet drew flow in through the branch, no longer reports merge coefficients under split names.

- `validation.junction.ideal_tee_score` scores `IdealMomentumTeeElement` against its pre-registration in #517, on the harness's imposed-flow networks. Against Bassett 2001's measured dividing data (role: selected) it over-predicts in every slice: K5 by +0.26, K6 by +0.23. The merge reproduces the scored `BM theory` variant to 0.8% wherever the branch stays below Mach 0.15. The harness's `MPCENetwork` gains a `_junction` hook, and `IdealMomentumNetwork` uses it.
- `IdealMomentumTeeElement` documents its range: subsonic supplying ports. A supplier that would really choke (Wang 2014, a = 2.44, high q and M3) is driven supersonic, and the solve stalls at Mach 1. It is not refused or flagged.

- C++ ideal momentum junction kernel: optional per-port receiver loss factor `IdealJunctionGeometry.loss_factor` (#526). It scales each receiver's vector Borda loss against the mixed supplied stream, `Pt_j = Pt_mix - f_j L_j`. `f = 1` (or empty) is the parameter-free closure, bit for bit, and any `f >= 0` keeps the mechanical dissipation non-negative. The element and GUI follow in later PRs.

- `IdealMomentumTeeElement(loss_factor=...)` (#526): one factor for all ports or one per port. It scales each receiver's loss against the mixed supplied stream; the suppliers' mixing loss is unchanged. The default 1.0 is the parameter-free closure, bit for bit. Any `f >= 0` keeps the mechanical dissipation non-negative; a split below `f ~ M^2/4` is refused by the energy check. The robustness harness (`validation/junction/random_robustness.py`) takes the factor as an optional fourth argument.

- GUI: the ideal momentum tee has a **Loss Factor** field (#526): it scales each receiving port's loss, with 1.0 the parameter-free model and a minimum of 0.25. The floor is from the robustness sweep: down to 0.25 convergence inside Mach 0.3 stays at or above 98.4%, while a nearly lossless split creates flow work and is refused. The help text calls it a rig-matching knob with no recommended value.

- **Tee port area vs the connected channel is reported.** Every multi-port junction now reports, per port, the connected element's own area and its ratio to the port area, plus `port_area_mismatch`: the number of ports that differ by more than 5%.
  - In the GUI, the tee Inspector shows an amber note naming those ports and both diameters.
  - Such a step is not modelled. The straight port takes the main pipe's area, so a narrower straight channel leaves an undeclared contraction there. The ideal momentum tee, which couples on static pressure, then creates total pressure across it, and its solve is rejected (#536).
  - Nothing is corrected; how to model the step is still open in #536.

- **GUI: the second-law verdicts and the network balance (#545).** Report only: a solve's success is unchanged.
  - **Status badge:** after a solve it says how many elements create total pressure, or that the balance does not close. It links to a "Network balance" window that lists what enters and leaves at each boundary, the heat terms, the mass and energy closure, and the failing elements.
  - **Inspector:** an amber "Second law" box appears on any node that owns a failing element. That includes a tee whose hidden port step fails, and a liner or array whose segments fail.
  - **Element results:** each carries its verdict (`second_law`) and `sigma_norm`, which is also selectable as a diagnostic.

- **Every converged solve reports a second-law verdict per element and a control-volume ledger (#545).** Report only: `__success__` is unchanged.
  - **`__second_law__`:** each element's internal entropy generation, `sigma = sum_out m s - sum_in m s - integral dQ/T`, with streams at their stagnation states and heat at the fluid's static temperature.
    - For an adiabatic two-port, `sigma_norm` is exactly `-ln(Pt_out/Pt_in)`.
    - Heating at Mach M costs total pressure and cooling may raise it; the heat term allows exactly that, bracketed over the element's Mach path.
    - The verdict is pass, fail or not_checked (with the reason). A crossflow chain is checked as one control volume; the legacy `TeeJunctionElement` is not checked.
  - **`__ledger__`:** mass and energy over the network's control volume (streams at the flow boundaries, heat given, heat leaving or withheld) and their closure.
  - **Measured on the test suite** (25,356 converged solves):
    - No element of a working network fails. The only failures are junction roots the existing junction check already rejects.
    - Re-creating the #541 area-change bug fails exactly that element.
    - The ledger closes to 1.1e-7, except where the solver's temperature clamp discarded energy, which it shows.
  - **C++:** `entropy_generation` and `heat_entropy_flux` (variable cp, exact at rest, Mach part by quadrature).

- A selectable C2 re-fit of the NASA-9 species thermo (#549): `set_thermo_model("nasa9_c2")` (C++ `set_thermo_model(ThermoModel::NASA9_C2)`). NASA-9 is C0 in cp at its interval joins but not C1: dcp/dT jumps 0.15-100 % at 1000 K and up to 48 % at 6000 K, and every Jacobian using it steps there. The re-fit bridges each join inside |T - T_join| < 0.05 T_join with a degree-7 cp/R that is C2 with NASA-9 at both window edges and whose h and s rises equal NASA-9's. Outside the windows it is identical to NASA-9; inside, cp moves by at most 4.8e-4, h by 1.1e-5 cp T and s by 3.1e-5 R. The bridges are generated by `thermo_data_generator/generate_c2_refit.py` into `include/thermo_c2_refit_data.h`. The default stays NASA-9.

- **Fanno-Rayleigh duct kernel (#551): friction and heat transfer together.** `fanno_rayleigh_channel_flow` has `fanno_channel_flow`'s contract, plus the exit stagnation temperature and a new `dG_dTt_exit` field. It marches the duct with Tt rising linearly from inlet to exit, exactly for a thermally perfect gas.
  - **Choking:** at the exit, at an interior throat (cooling after a sonic point), or at the inlet.
  - **Also new:** `fanno_rayleigh_march` (one duct at a given mass flux) and `fanno_rayleigh_choked_mass_flux`.
  - **Checked against:** the adiabatic kernel at zero heat (<= 1e-7 in G), and exact Rayleigh flow at zero friction (<= 1e-8 in exit Mach), including its thermal choke.
  - **Not yet wired into `ChannelElement`;** a follow-up does that.

- **A wall-heated compressible channel flags the heating its momentum leaves out (#551).**
  - **Why:** its kernel is adiabatic Fanno flow at the inlet Tt, and the wall heat only raises the downstream node's Tt.
  - **New diagnostics:** `rayleigh_omitted_frac` and `rayleigh_dPt_omitted` give the first-order Rayleigh estimate (gamma M^2 / 2) ln(Tt_exit / Tt_in) of the missing total-pressure change: a loss when heated, a gain when cooled. `rayleigh_omitted` is set above 0.1% of the inlet Pt.
  - **GUI:** the Inspector shows "Heating not in the momentum" with the estimate.
  - **Size:** in the stress hotcold pairs it is +4 to +9% of Pt on the heated side and -5 to -24% on the cooled side near choke.
  - **The solve itself is unchanged.**

- GUI: the solver settings take the incompressible validity limit ("Incompressible Mach Limit"; blank = 0.2), and an ejector's Inspector shows Sanger's four jet-pump loss coefficients when the network is incompressible. Blank nozzle / suction = the compressible nozzle efficiencies, K = 1/eta - 1. Its result block shows P_c* only where it exists, plus the jet pump's head ratio, efficiency and jet Mach (#588).

- `MomentumChamberNode` takes a `regime` (#588), default `DEFAULT_REGIME` like the elements, and the GUI gives every chamber the network's global regime. `"incompressible"` is the low-Mach limit of the chamber closure: `T_static = T0` and `Pt = P + (m/A)^2 / (2 rho)`, with no sonic limit. The merge face follows, linear at constant density. Both match the compressible forms to O(Ma^2) (0.25% of the head at Ma 0.2), so an incompressible network no longer hides a sonic limit in its chambers. A script that relied on the default chamber being compressible must now ask for `regime="compressible"`.
- The incompressible validity report (#592) covers these chambers too: a chamber on its incompressible closure above the Mach limit is listed by its node id, and the GUI attributes a junction's inserted port chamber to the junction.

- `EjectorElement` takes a `regime` (#588), default `DEFAULT_REGIME` like every element, and the GUI passes the network's global regime. `"incompressible"` is the ejector's low-Mach limit: Sanger's one-dimensional jet pump (NASA TN D-4445, 1968), with each stream at its own constant density, no choking and no critical regime.
  - Optional loss coefficients `k_nozzle`, `k_suction`, `k_throat` and `k_diffuser`. The nozzle and suction defaults are the compressible closure's own nozzle efficiencies (K = 1/eta - 1), so both regimes agree at low Mach out of the box.
  - It reproduces Sanger's Eq. 11 to machine precision.
  - With his calibrated losses, his theory scores 9.3% on his own measurements (TN D-4759, Table II).
  - The compressible ejector reduces to it as Ma^2.
  - A script that relied on the default ejector being compressible must now ask for `regime="compressible"`.

- Every converged solve reports each element running an incompressible closure above one Mach limit (`NetworkSolver.incompressible_mach_limit`, default 0.2), judged at its own highest internal Mach: an orifice's jet, a channel's higher end. Report only: `__incompressible_validity__` and `__incompressible_warnings__` in the result, a line on the GUI's status badge and a block in the Inspector. At Ma 0.2 the measured error is 4.0% on a channel's drop and 2.9% on an orifice's; the limit is a rule of thumb, not a guarantee (#588).

- **Continuation in the boundary conditions, after the Newton homotopy (#642).** `NetworkSolver.bc_ramp_retry` (opt-in; the GUI, `NetworkRunner` and the stress harness turn it on) gives a cold solve that the homotopy did not rescue one more try. Every boundary is ramped from its gentle end, where every pressure boundary sits at the mean of their Pt and every flow and heat boundary at zero, up to its own value, each step warm-started from the last root.
  - On the stress random draws the homotopy left unsolved, it rescues 3 of 21 (in 13-46 s).
  - Raise the max solve time for networks that need it.

- `discharge_pipe_zeta(correlation, hole, beta)`: Idelchik (1966) diagrams 4-13..4-16, the permanent loss of a discharge hole as an orifice in a straight pipe, on the pipe velocity head. Its zeta' and tau come from the wall diagrams' own tables.
- `discharge_pipe_cd_ratio(correlation, hole, beta)`: the factor that carries a plenum-measured Cd into the pipe through ISO 5167-2's recovery relation. It tends to 1 as beta tends to 0.
- `OrificeElement` diagnostics report `Cd_pipe_ratio` along the flow.

- An `EffusionLiner` is a film provider: a `ThermalWall` downstream on the gas path can name it (`film="<liner id>"`). Its panels are one group of Baldauf rows each, at their own blowing ratio, density ratio and coolant temperature, layered row by row (`cb.film_span_average_baldauf_2002_groups`, new; `FilmSpanAverageGroups`). Every panel's flow and station temperature reach the wall's Jacobian. `EffusionLiner(film_direction=...)` orders the panels along the gas path. In the GUI the liner node has the "f" film handle and a gas-direction select, and is listed in a wall's Film dropdown (#695).

- Slot film-cooling effectiveness, Marek & Tacina (1975) NASA TN D-7958: `slot_film_effectiveness_marek_tacina_1975` and its `_and_derivatives` (eta, deta/dM) form. Juhasz-Marek mixing with the mixing coefficient from the free-stream turbulence (the report's C_m = 0.25 Tu/100 with Tu in percent; `Tu` is passed as a fraction, 0.14 = 14%, like every turbulence input, and above 1 raises) and a user `C_m_multiplier`. The first slot (liner cooling ring) film kernel; the hole films do not apply to a slot (#714).

- `SlotFilmElement`: a 2-D film slot (a combustor liner's cooling ring) as a film provider for a downstream wall. Film: Marek & Tacina (1975), C_m from the free-stream turbulence `turbulence_intensity`, a fraction as on every element (default 0.14 = 14%, the D-7958 case Eq. 5 fits best), with a user `C_m_multiplier`. Plus `slot_film_span_average_marek_tacina_1975`, the span average a wall sees (#714).

- GUI: a **Film Slot** node (`slot_film`), the `SlotFilmElement`, in the palette and the Inspector. Defaults are the D-7958 rig (0.63 cm slot, 37 cm span, Tu 0.14 = 14%). A thermal wall downstream names it as its film, like film holes (#714).

- `ImpingementArray.summarize` reports starved rows as unresolved: `rows_unresolved` and `n_rows_unresolved` mark a row whose plate total-pressure drop is below 1% of the gap's dynamic head, where the crossflow chain cannot resolve the jet flow (its sign included) or the Nu and eta_r that follow from it (#729). On NASA CR-3630 Table D.2 it marks exactly the 8 starved rows.

- **Area change closures exposed (#731):** `sudden_expansion_zeta` (returns `(zeta, dzeta_dRe)`), `sudden_contraction_zeta` (returns `(zeta, dzeta_dRe)`), `conical_diffuser_zeta`, `conical_nozzle_zeta`, `compressible_borda_carnot`, `area_change_loss`, `area_change_flow` and `area_change_element_flow`.
- **`viscosity_dX` (C++),** the analytic composition derivative of the Wilke mixture viscosity.

- **The effusion plate reports what its gas-side augmentation multiplies (#476a).** Beside `h_gas_unblown` the diagnostics now give the chamber baseline's model (e.g. `SmoothModel(gnielinski)`, fully developed), Re, Pr and Nu, and the Inspector shows them under *Gas h Augmentation*. On Andrews' rig the required augmentation moves ×1.75 between two baselines, so an augmentation taken from the literature only carries over against the same baseline.

- `EjectorElement(anchor_closure=...)` selects the pair of operating anchors. `"huang1999"`, the new default, is Huang 1999's printed P_c* chain (φ_m mixing, normal shock, isentropic diffuser) with a sudden-expansion dead head. On the measured air rigs it halves the error of the whole entrainment curve (12-23% RMS against 29-48%), and no stress root it returns breaks the second law.
  - `"kracik_dvorak"`, the lossless pair of 0.6.0, stays selectable.
  - The Python reference and analytic C++ Jacobians are `huang_critical_back_pressure`, `sudden_expansion_dead_head`, `ejector_huang_critical_back_pressure_and_jacobian`, `ejector_sudden_expansion_dead_head_and_jacobian` and `ejector_huang_mixing_coefficient`.
  - The GUI Inspector shows it as "Back-Pressure Anchors".
- A failed solve whose ejector outlet lies above the dead-head back pressure now says so in the solver message: a real ejector backflows there, which the model does not cover.

- GUI: a thermal wall's box can be dragged off whatever its auto placement put it under. The offset is saved with the wall and follows it when its nodes move; double-click the box to reset it. Moving it keeps the solved results. Wall boxes now draw above unselected nodes, and a selected node still comes to the front.

- `velocity_of_approach_factor(beta)`: ISO 5167-1's E = 1 / sqrt(1 - beta^4), with 1 outside 0 < beta < 1, the network's no-correction convention. It is now the only definition: the orifice kernels' six inline copies call it.
- `OrificeElement` diagnostics report `beta_flow` and `E`, the velocity of approach along the solved flow. A reversed orifice reports its to-node pipe's; `beta` stays the declared one.
- GUI orifice Inspector:
  - The **Pressure Taps (reported)** select (Flange, Corner, D and D/2) for the ISO 5167 correlations. It picks the differential reported as dp at taps and never moves the flow.
  - The reverse-flow approach pipe with its beta and E.
  - The solved approach-pipe static pressure p1, the tap differential, the permanent loss and its ratio, and a warning when the tap pressure falls outside ISO's permanent-loss range.
  - beta and E now come from the solver's resolve pass rather than being recomputed in the browser.

- **Element diagnostics report where each geometry value came from (#559).**
  - Every element that inherits or derives its flow area now surfaces `area_source` in its `diagnostics()` -- `"user"`, an `"inherited from ..."` / `"channel '<id>'"` string, or `"default"` -- and an area change reports `F0_source` / `F1_source` per face. `ChannelElement` and `TeeJunctionElement` previously discarded this provenance; `MomentumChamberNode`, `CombustorNode`, `PressureLossElement` and `AreaChangeElement` computed it but never exposed it.
  - This is the authoritative provenance the GUI can read instead of re-deriving inheritance itself -- groundwork for single-source-of-truth geometry inheritance.

- **Orifice and junction-port geometry provenance (#559).**
  - An `OrificeElement` now reports `area_source` (`"user"` or `"default"`) plus its bore `diameter` and `beta` in `diagnostics()`; the fallback bore is the named constant `DEFAULT_ORIFICE_BORE` (0.08 m) instead of a bare literal.
  - A `MultiPortChamberBase` reports `port_<i>_area_source` (`"user"` or `"inherited from <id>"`) per port, set idempotently across resolve passes.
  - Completes the geometry-provenance exposure started in B1a (#560).

- **The junction closure is checked against its own paper, Mynard & Valen-Sendstad 2015.**
  - **Before:** the C++ port was held to the Python one, and nothing held the Python to the paper.
  - **Fidelity:** the kernel in Mynard's configuration draws his model's curves on all 31 panels of Figs 4 and 6-11, within 4.2 px of the 150-dpi scans and 2.7 px on 29 of them (<= 0.032 in K). Fig 4's eta_j = 0 curve matches too. The curves were digitised by colour, reproducibly, with `validation/junction/digitise_mynard2015.py`.
  - **His CFD (laminar, Re 350-2400)** is reported against both configurations by `python -m validation.junction.mynard_fidelity`, and is not tuned on. Over all 158 samples: MAE 0.212 for his model, 0.219 for production.
  - **New argument:** `_mynard2010.junction_loss_coefficient(dividing_streamline_recovery=...)`. The recovery is a repo addition with no switch until now, so `eta_scale=1.0` alone never gave his model. The default is unchanged.
  - **Bassett 2001's Table 2 transcription** is now pinned against his own calculated curves.

- **Pressure taps for the metering orifice correlations.** `OrificeGeometry.tap` (`OrificeTap.Corner`, `D_D2` or `Flange`, the default) selects the ISO 5167-2:2003 tap positions for Reader-Harris/Gallagher and Stolz.
  - **What it changes:** the tap is a measurement convention. It changes which pressure difference C refers to, not what the plate does to the flow, so the network element evaluates fixed taps and its flow never depends on the tap (D-D/2 since #645; `OrificeElement(tap=...)` picks only the reported differential).
  - **New functions:** `metering_cd_and_dRe(correlation, geom, Re_D)` returns (Cd, dCd/dRe_D). `metering_cd_in_range(...)` flags ISO 5167-2's d, D, beta and Reynolds limits, and the orifice's `Cd_in_range` diagnostic now uses it.
  - **Stolz's tap default:** Stolz through `make_correlation` / `solve_orifice_mdot` now honours `geom.tap`. Before, it was always corner taps; at the flange default, Cd moves by about +0.0005 at beta 0.5 in a 100 mm pipe.

- **`/resolve_preview` endpoint: geometry inheritance without solving (#559).**
  - A new backend endpoint (backed by `gui.backend.resolve_preview`) runs the authoritative resolve pass without solving and returns, per drawn element/node, each geometry field's resolved value and provenance (`"user"` / `"inherited from ..."` / `"default"`), plus whether it is unresolved -- the network is not yet solvable, and why.
  - This is the single source of truth the GUI will read for inherited placeholders and to flag required values before Solve, instead of re-deriving inheritance in its own graph walk (B2 of the geometry-inheritance roadmap).

- A strict solve (`COMBAERO_SOLVER_STRICT=1` or `raise_internal_errors = True`; the test suite runs strict) now checks every residual row's declared units against its Jacobian at the converged solution and raises `RowKindMismatchError` on a misfiled row. `NetworkSolver.check_row_kinds(x)` returns the per-row report. A row's units are read from an exact unit coefficient (+-1 within 1e-3) on a pressure or a flow unknown. With the #614 fix reverted, 51 existing tests now fail naming the element; before, the suite passed with the bug in place. Lenient (production) solves are unchanged and pay nothing.

- **Second-law report: each ThermalWall as its own control volume (#551).** `__second_law_walls__` gives each wall's entropy generation `sigma_wall`. It is the entropy the heat carries into one stream minus what it took from the other, each integrated along that stream's own temperature path. It goes negative when an uncapped wall drives heat from the stream that is colder on average to the hotter one. The #551 seed-23 network "converged" with its 1600 K stream cooled to 298 K against a 400 K cold inlet; it now reports `fail`. Report only; `__success__` is unchanged. A check at the temperatures the wall model sees would be vacuous (`Q = U (T_A - T_B)` makes it non-negative by construction), so it is not used.

- **`fanno_channel_flow(..., dY_directions, with_derivatives)`** (#489).
  - **`dY_directions`:** mass-fraction directions v; `dG_ddir` returns dG/ds along each (Y -> Y + s v). A unit vector e_k gives dG/dY_k with the other fractions fixed. Exact to 1.5e-7 against Richardson differences on the unchoked, choked and near-choke branches, and to <= 3e-3 of the column at 3 Pa of drive on 2 bar. Where a central step would make a fraction negative (a species absent from the node but carried by an inflow), the difference is one-sided.
  - **`with_derivatives=False`:** returns G, choked and the Mach numbers only, for residual-only evaluations.

- **The compressible `ChannelElement`'s Jacobian carries the feeding node's composition (#489).** Where the plenum feeding a duct mixes gases, the composition moves with the inflows, and that column was missing: the inflow columns of the global Jacobian were 1.6-2.9% off. They now agree with differences to ~1e-8.
  - **Directional:** the solver passes an orthonormal basis of the directions the node's composition moves in. That basis comes from an SVD of its relay dY/dx and is one vector for two mixing gases. The kernel differentiates along it, and the gradient is projected back onto per-species keys, which is exact for every relay column. A per-species gradient made the column cost more than it saved: 2.2x the solve time eagerly, still +11% once lazy.
  - **Lazy:** an element's `residuals` may now return a `LazyJacobian`, a mapping computed on first read. The compressible channel does, so its implicit Fanno derivatives and the column are computed only when the solver asks for J, and residual-only evaluations skip them.
  - **Measured** on 36 H2/CH4/CO2 mixing networks: 431 -> 401 evaluations and 2.3 -> 1.8 s against main, all converged as before. Solutions agree with main's to a median 4e-12, max 5.5e-8.

- **`transport_and_dT`, `omega22_and_derivative`, `dcp_R_dT`:** mixture viscosity and conductivity with analytic d/dT, the collision integral with its slope, and per-species d(cp/R)/dT (#485). At the NASA polynomials' 1000 K range join, cp's slope jumps, so dk/dT is one-sided there.

- **`fanno_channel_flow`:** the mass flux a duct passes from a stagnation state against a back pressure, with dG/dPt0, dG/dTt0 and dG/dP_target. It saturates at the choked flux and is never infeasible.
  - **Method:** it solves for the exit Mach, because P_exit is linear in G at fixed M so G(M_exit) is explicit, and differentiates implicitly. That stays well conditioned up to the choke and C1 through it.
  - **Below 1 Pa of drive** the flux is blended C1 to a finite slope.
  - **Accuracy:** it inverts `fanno_duct` to 1e-11.
  - **Cost:** about 4 ms per unchoked call, 6 ms choked.

- **Fanno flow integrated in Mach number** (`fanno_mach.h`: `fanno_duct`, `fanno_choked_mass_flux`, `fanno_length_between`, `fanno_inlet_mach`, `fanno_sonic_mass_flux`, `fanno_state_at_mach`, `fanno_dx_dmach`). It is the kernel for the compressible-channel redesign.
  - **Method:** at fixed mass flux the station state is algebraic in M, so the length is a smooth integral that is exactly 0 at sonic. L* and the choked flux are exact, with no M = 0.999 cutoff, gradient floor or step control.
  - **Accuracy:** it agrees with the 20000-step x-march to 1e-13 to 1e-9 at constant f and to ~1e-8 with Haaland. Its choked flow matches the march's choke bisection.
  - **Cost:** 1.5 ms per duct, 2.8 ms per choked-flux root.
  - **Why not the march:** near choke the x-march stalls at 1e-4 to 1e-5 whatever its step count.

- **Effusion liner: a duct-fed effusion wall (#471).** Python `EffusionLiner` and the GUI node "Effusion Liner": coolant in, coolant out and discharge ports.
  - **Structure:** N stations along the backside duct, each a Bassett bleed (`CrossflowSegmentElement`) feeding a duct-fed panel into the gas.
  - **Crossflow Cd:** panels take McGreehan-Schotsch's supply-side U1/Vi from the duct's mean velocity at their station and the static-referenced ideal jet velocity (C++ `crossflow_velocity_ratio`). The derivative chain runs into the neighbour segments.
  - **Physics as solved:** the backside static pressure rises along the duct, the downstream panels pass more, and the upstream panels (faster crossflow) have the lower Cd.
  - **Flags:** U1/Vi beyond the Rohde-scored bounds (0.2, 0.35); psi beyond Bassett's measured 3; the provisional plenum-fed coolant-side heat transfer.
  - **GUI default:** Andrews plate C holes over a 30 mm duct. With a 200 Pa duct drop against ~3.8 kPa hole drive, U1/Vi peaks at 0.27.
- **`EffusionPlateElement(crossflow_segments=..., crossflow_area=...)`:** duct-fed panels. `OrificeElement` discharge-hole Cds take `U1_over_Vi`, and their Jacobian carries every Cd dependence by name.

- **Bleed through a duct wall: `CrossflowSegmentElement` with Bassett stations (#471).** This generalises the impingement crossflow segment.
  - **Stations:** a channel segment whose end nodes are side-stream stations, merge or bleed. The C++ (f, J) `station_half_drop(m_a, m_b, P, T, X, A, kappa)` carries half of each station at its node's density, and `channel_entry_drop` covers a reservoir entry.
  - **Merge:** `kappa = 0` (normal injection) is Florschuetz's term.
  - **Bleed:** `kappa = 0.75` is exactly Bassett, Winterbone & Pearson's (2001) separating straight-run K2/K5 (Eq. 15) for every split. Scored on Bassett's measured K5 (Fig. 7c, 43 points): MAE 0.039 in K, 0.021 at q >= 0.6.
  - **Inheritance:** `ImpingementCrossflowElement` is now a `kappa = 0` subclass.

- **`EffusionPlateElement` owns its wall; the discharge node decides the gas side (#471).**
  - **The wall.** The plate is the wall the coolant passes through. Its diagnostics now carry the wall solution: `T_wall_hot`, `T_wall_cold`, `q_wall`, `Q_wall`, `eta_overall`.
    - **Coolant side:** Andrews 86-GT-225.
    - **Conduction:** `wall_thickness / wall_conductivity`, through C++'s wall coupling.
  - **Gas side, by discharge node:**
    - **A momentum chamber (flow):** its own surface correlation gives the unblown coefficient, and its velocity gives the blowing ratio. The gas temperature is the approaching main stream's, not the mixed chamber state, which the effused coolant dilutes. `gas_augmentation` is the caller's input, default 1.
    - **A plenum (state, no flow):** the imposed `gas_heat_flux`, default 0 (adiabatic).
  - **Network energy:** unchanged. The heat returns with the coolant into the same node.
  - **Film:** off by default. Baldauf + Sellers (`gas_film="baldauf_sellers"`) is selectable and flagged outside its envelope. Its panel average moved into C++ (`effusion_panel_film_effectiveness`), shared with the Andrews validation runner, which gives bit-identical results.
  - **Coolant crossflow:** a channel-fed supply is flagged (`coolant_crossflow_ignored`). This 2-port plate is plenum-fed.
  - **One wall core:** `overall_effectiveness` now uses the same C++ wall core, with zero conduction resistance. Its results are unchanged.
  - **GUI:** new "Effusion Plate" node, defaulting to Andrews' plate C, the geometry the closure is scored on. Connect it to a chamber's side handle.

- **The merge chamber, the jet momentum and the impingement crossflow term are compressible (#471).** The incompressible first version used one density for the main face and the outlet, `Pt - P = ½ρu²`, and a jet density at the supply total temperature. All three evaluations now come from C++:
  - **Exact impulse.** `chamber_merge_face_state` solves the face from the exact impulse balance (the ideal-gas quadratic on its subsonic root, flagged when no such root exists). It takes the face Pt from the chamber's own entropy-based closure.
  - **Jet momentum.** `jet_impulse` uses the isentropic jet velocity, or the sonic momentum plus pressure thrust when choked.
  - **Crossflow segment.** The crossflow term puts each station at its own density (now `station_half_drop`).
  - **Jacobians.** The density derivatives come from `density_and_jacobians`, so no density derivative is computed in Python.
  - **Validation shift.** The impingement chain's Gc/Gj vs Eq. 8 (an incompressible model) moved from max 3.7% to 4.1% over a 3% pressure drop. Fig. 6 is unchanged (+4.04% / 6.89%).

- **Merge chamber: a `MomentumChamberNode` with a declared main inlet accepts side streams (#471).** It still has one outlet and is closed by the axial momentum balance over a constant-area control volume, with no loss coefficient.
  - **Momentum:** each side stream brings `m u_jet cos(theta)` of axial momentum and discharges at the chamber's static pressure. The main-inlet element sees the main-face state (`cb.chamber_merge_face_state`).
  - **Validation:** the solved network satisfies the control-volume balance, normal injection costs the main stream its mixing loss, and inclined jets pump it. The global Jacobian matches finite differences to 1e-4.
  - **Injection angle:** orifice-type elements gain `injection_angle_deg` (90 by default). `EffusionPlateElement` uses its hole angle.
  - **Area:** inherited from the main inlet; a mismatch is refused.
  - **Splits:** still refused.
  - **GUI:** the chamber node gains a side-stream handle.
  - **Solver:** merge networks join the LM fallback, restarted from the cold guess.

- **Impingement array, channel exit, as one element (#465): Python `ImpingementArray` and the GUI node "Impingement Array".** It models exactly one configuration: a jet plate fed from the upstream node, impinging on a target, with the spent air leaving down the gap to the downstream node (Florschuetz et al. 1981).
  - **Expansion:** the array becomes the validated row chain (`{id}__p{i}`, `{id}__x{i}`, `{id}__c{i}`).
  - **Thermal wall:** one wall to it is spread over the rows, each row on its own target footprint.
  - **Results:** reported back under the array's own id:
    - total flow and supply-to-exit dP;
    - Re_j range, max Gc/Gj, mean Nu/h, jet-flow non-uniformity;
    - out-of-range flags;
    - per-row lists, shown as a table in the Inspector.
  - **Wall result:** the edge reports the total Q and the hottest row's wall profile.
  - **Default:** a Florschuetz test geometry (0.1 in holes, (5, 4, 2), 10 rows) that solves inside the correlation's range.
  - **Other configurations:** a bypass crossflow (#467) and spent air leaving through the target (#468) will be separate elements.
- **GUI orifice offers `Lichtarowicz`** (long orifice, l/d 2-10).

- **Generic C++ core for jet-plate rows (#465).** The physics behind
  `ImpingementPlateElement` and `ImpingementCrossflowElement` moved from
  Python into C++ (f, J) functions that later configurations reuse:
  - `jet_row_heat_transfer(set, JetRowGeometry, m_jet, m_crossflow, mu, k, Pr)`
    returns h with analytic `dh_dm_jet` / `dh_dm_crossflow`;
  - `side_stream_momentum_drop` (since replaced by `station_half_drop`, #471) gives the centred
    momentum term for side streams joining a channel.

  The validation numbers are unchanged to the digit.

- **Impingement-cooled walls as a network (#465).** Two new elements, one
  pair per spanwise row, model a plenum-fed jet plate with its target and
  crossflow channel.
  - **`ImpingementPlateElement`:** an orifice of `n_holes` holes whose
    flow is the row's jet flow, plus Florschuetz et al. (1981) heat
    transfer on the target. **Gc/Gj comes from the network:** the crossflow
    approaching the row is every other inflow to the plate's crossflow node,
    not the uniform-supply closed form.
  - **`ImpingementCrossflowElement`:** a channel carrying the momentum the
    jets cost, `(m_b|m_b| - m_a|m_a|)/(rho A^2)` per merge, split between
    adjacent segments.
  - **Result:** over the 27 Fig. 6 geometries the network reproduces
    Florschuetz's own flow model (Eq. 8) to 3.7% max in Gc/Gj. Without the
    momentum term the supply stays uniform and misses it by 2.2x.
  - **Fig. 6 Nu/Nu1 through the chain:** +4.0% bias and 6.9% MAE, against
    the correlation's +3.7% / 6.8%.
  - **Solver opt-in:** an element defining `network_flow_inputs()` receives
    the neighbours' flows as `flows=` on `residuals`, `htc_and_T` and
    `diagnostics`. An htc that reads them returns `dh_dsources`, which the
    wall relay carries into the Jacobian.
  - **C++:** `JetArrayImpingementResult.dNu_dGc_Gj`, analytic.
- **`'Lichtarowicz'` hole Cd on `OrificeElement` and `EffusionPlateElement`
  (#465).** Lichtarowicz, Duggins and Markland's (1965) long-orifice
  correlation (l/d 2-10, Re 10-2e4) already existed in C++ but no element
  could select it. Below l/d = 1.5 it is refused at set-up, naming the
  element.
- **`discharge_cd_in_range` and the `Cd_in_range` orifice diagnostic (#465).**
  Reports whether a discharge-hole Cd is inside its source's own Re and l/d
  range. It is a flag, never a selector.

- **Chyu (1990) pin-fin sets, and inline friction (#335).**
  - `chyu_1990_nu(arrangement, fillet)`: Table 2, straight or pin-endwall-fillet pins, inline or staggered (pin surface).
  - `chyu_1990_friction(arrangement, fillet)`: Fig. 6, digitised and fitted (provenance Fitted). The points are committed in `validation/cooling/data/chyu1990/`.
  - `chyu_1990_fillet_over_straight`: a fillet modifier (staggered 0.78-0.92).
  - **Inline pin-fin arrays now solve.** `PinFinModel` defaults inline friction to Chyu's set, the only inline source, measured at one geometry. A staggered friction set on an inline array is refused.
  - **New friction basis:** `PinFrictionBasis.PerRowHalfHeadVmax` for Chyu's `f = 2 dP/(rho Vmax^2 N)`, 4x the canonical value.
  - **Cross-lab disagreement:** Chyu's staggered friction is 16-25% below Metzger's at the same geometry. Reported, not reconciled.

- **Pin-fin channels in the network and GUI (#335).** `PinFinModel` on a
  `ChannelElement`'s `ConvectiveSurface` uses the pin-fin sets.

  - **Velocity:** Vmax comes exactly from the array geometry.
  - **Heat transfer:** `h` is referenced to the endwall's base area and
    includes the fin efficiency of pins fed to mid-height (`k_pin`, default
    20 W/(m K)).
  - **Drop:** the correlation's own per-row `dP = 2 rho Vmax^2 N f(Re_D)`,
    with an analytic Jacobian in mass flow and upstream T and P (through rho
    and mu).
  - **Inline arrays require a user friction set.** None is available, and
    staggered friction is not substituted.
  - **GUI:** the "Pin-Fin Array" surface type is back, and networks saved
    before 0.7.0 load under their old field names (`pin_diameter`,
    `channel_height`, `is_staggered`).
  - **Known gaps:** `T_aw` is borrowed from the smooth correlation, and
    `dh_dT` is a scaled smooth-correlation proxy. Both are documented.

- **Pin-fin arrays: provenanced sets, exact basis converters and a fin
  model (#335).** `pin_fin_correlation.h` carries each source in its own form
  and converts it exactly, from geometry alone, to one canonical basis:
  `Re_D` on the pin diameter and the velocity at the minimum flow area,
  `Nu_D = h D/k`, and `f = dP/(2 rho Vmax^2 N)`. Heat transfer and friction
  are separate sets.

  | set | form | box |
  |---|---|---|
  | `metzger_1986_staggered_nu` (default) | `0.135 Re^0.69 (X/D)^-0.34` | H/D <= 3, 2 <= S/D <= 4, 1.5 <= X/D <= 5 |
  | `metzger_1982_staggered_friction` (default) | two branches, C1-blended at Re 1e4 | 0.5 <= H/D <= 6, 2 <= S/D <= 4 |
  | `vanfossen_1982_staggered_nu` | D' basis | H/D 0.5-2 |
  | `damerow_1972_staggered_friction` | per (N-1) rows | S/D 4.24-7.07 |
  | `chyu_1998_nu` | staggered/inline x pin/endwall/total | S/D = X/D = 2.5, H/D 1 |

  - **Converters:** D'/D, A'/A_min, D_h/D and A_min/A_frontal (the smaller of
    the transverse and staggered diagonal gaps), plus per-wall area
    fractions.
  - **Fin model:** `pin_fin_array_efficiency`, with an analytic d/dh.
  - **Data-backed ratio modifier:** `chyu_1998_inline_over_staggered`,
    0.76-1.0 over Re 5k-30k.
  - **Inline friction:** no source is in hand, so none ships, and the
    modifier says so (`has_f = False`).
  - **Behaviour:** outside its box a set is flagged `extrapolated`, never
    refused. Evaluators never throw, and derivatives are analytic.
  - **What the code removed in #332 got wrong, confirmed against the
    sources:**
    - dP was a factor of 4 low (`rho V^2/2` against Metzger's
      `2 rho V^2`);
    - the friction exponent was -0.25;
    - a Pr^0.4 term and three spacing/length terms were invented;
    - friction was credited to Simoneau & VanFossen.
  - **Sources:** see `validation/cooling/extractions/pin_fin_sources.md`.
  - **Not yet included:** the network element (`PinFinModel`) and validation
    scoring follow.

- **Taslim & Spring (1987), the first shipped ratio-form rib sets (#444).**
  `taslim_spring_1987(aspect_ratio_taslim, e_D)` returns one `RibRatioSet`
  per configuration the paper tested with two ribbed walls (seven: AR 0.5 /
  1.0 / 3.5, transverse ribs, p/e 10), in the form the paper states:
  `Nu/Nu_DB = C (Re/1e4)^-0.2` (their `Nu ~ Re^0.6`) and a Re-independent
  passage-average Fanning `f`. `C` is read from figure 9's markers and `f`
  is figure 11's mean, so provenance is `Fitted`; an untested configuration
  throws rather than interpolating. Taslim's aspect ratio is height/width,
  so the set's `W/H` is `1/AR`.

  - **New field `RibRatioSet.Re_floor_f`:** friction's own lower bound
    (0 = use `Re_floor`). Taslim measured friction to lower Re than heat
    transfer; with one shared floor, f drifted off the paper's constant
    inside its own measured range.
  - **Scored on figures 4/5**, the paper's own Nu-vs-Re data, so this is
    fidelity: 43 points, all in range, pooled MAE 6.1%, per configuration
    2.6-10.7%. The two largest (AR 0.5 e/D 0.125 at +10.7%, AR 1.0 e/D
    0.167 at -9.6%) are offsets between figure 9 and figures 4/5 in the
    paper itself. Dropping the `Re^-0.2` raises the MAE to 45%, so the
    score responds to the form.
  - **Friction is not re-scored:** figure 11's AR 1.0 series stay with
    `han_1988_orthogonal` as an independent-lab check (#443); the ratio
    sets' f comes from that same figure, so scoring it there would be
    circular.
  - Validation harness: a ratio runner beside the rib/jet/orifice runners,
    and `y_scale` in series metadata for figures plotted as `Nu x 1e-3`.

- **Ratio-form rib correlation sets: `Nu/Nu0` and `f/f0` beside `R`/`G`
  (#444).** `RibRatioSet` + `evaluate_rib_ratio` carry a correlation as a
  multiplier on the smooth-duct baseline its source fitted to, for sources
  that publish ratios -- without forcing their data through the
  law-of-the-wall similarity `R`/`G` assumes. The type only: which
  correlations to ship comes next, with Taslim & Spring and NASA CR-4396 as
  acceptance data.

  - **In range** `Nu = r * Nu0_source`: the paper, whatever baseline a
    caller prefers (re-referencing exactly cancels).
  - **Below the fitted Re floor** a C1 blend in `ln Re` hands over to
    `k * Nu0_ext`, `k` matched at the floor; `Nu0_ext` is smooth-pipe
    Gnielinski by default, so Nu keeps its laminar 3.66 limit as Re -> 0
    instead of collapsing with Dittus-Boelter. Friction holds its ratio at
    the floor.
  - **Robust for user coefficients too:** never divides by `Nu0`; `Re`
    enters as `sqrt(Re^2 + 1)`, so Nu and f are even and their analytic
    derivatives odd and finite through `Re = 0`; mistakes are rejected by
    `validate_rib_ratio_set` / `validate_rib_ratio_options`, and the
    evaluator never throws.
  - The handover uses smooth-pipe Gnielinski's exact derivative (see the
    #446 entry). Its blend anchors (Re 1000, 2300) are now named constants
    in `heat_transfer.h`.

- **Rib shape, and Han, Zhang & Lee (1991)'s nine configurations -- closes
  #434.** `RibCorrelationSet` gains `shape` (`RibShape`: Transverse,
  Parallel, Crossed, V, Lambda; `Unspecified` binds nothing), and the rib
  runner refuses a series of another shape for the SHAPE, not incidentally
  through rib angle. `han_zhang_lee_1991(shape, alpha_deg)` returns the
  paper's Table 2 row -- verified on four channels, including a reviewer
  reading the original -- for Transverse at 90 deg and Parallel, Crossed, V
  and Lambda at 45/60 deg, one rig (square, `e/D` 0.0625, `P/e` 10);
  anything else raises. Sets can now carry a source's OWN printed `G_bar`
  (`C_Gbar`, `RibResult.G_bar`/`has_G_bar`); the runner uses it instead of
  Han's 1.2 wherever it exists. V and Lambda sets are not `symmetric`.

  | scored, as fidelity | n | MAE |
  |---|---|---|
  | figure 4.51 crossed / V / Lambda, `G` | 35 | 0.9-2.9% per class |
  | figure 4.51 crossed / V / Lambda, printed `G_bar` | 26 | 1.9-2.4% |
  | figure 4.53 60 deg V, `Nu_ratio` vs `f_ratio` (accuracy) | 5 | 2.5% |

  66 of #434's 71 points score; figure 4.54's 45 deg BROKEN V ribs (4)
  stay refused with the reason recorded -- no source on disk covers broken
  ribs. The four earlier sets declare their shape and are otherwise
  unchanged.

- **Independent-lab rib friction is scored, and it closes #403 item 2
  (Lau vs Han) without picking a winner.** `runner.py` gains a `Re_Dh`
  path for passage-average friction: Taslim & Spring (1987)'s three AR 1.0
  two-side-turbulated series -- Han's own square-channel configuration from
  Northeastern / GE -- score against `han_1988_orthogonal` in `f` directly.
  The set's four-sided `f` is converted with Han's own `fbar` decomposition,
  so nothing on the measured side passes through `R`.

  | `e/D` | Han `fbar` vs measured |
  |---|---|
  | 0.083 | -14.9% |
  | 0.167 | -24.9% |
  | 0.250 | -49.1% |

  All out of domain (every `e/D` above Han's 0.047-0.078). On Han's basis
  Taslim's `R` is 2.77 against Han's 3.2 and Lau's ~3.6: the three labs
  BRACKET Han by about +/-13%, so there is no outlier to arbitrate. That
  spread is now the recorded cross-source accuracy of the Han family's
  friction. A test falsifies the convention trap: read as the four-sided
  `f`, the comparison flips to ~+60%.

- **Figure 4.53's 90 deg continuous rib is scored on the performance curve
  (#435).** `runner.py` gains an `f_ratio` path: `Nu_ratio` against
  `f_ratio` with Re eliminated, so Re is recovered by bisecting the set's
  own chain until its `f_ratio` lands on the abscissa -- the same move the
  `e+` path makes. Against `han_1988_orthogonal`: 5 points, MAE 7.9%, bias
  -2.1%, labelled ACCURACY (Han and Zhang 1992 is a different study).

  **Two conventions, and only one of them matters.** `f_ratio` is
  CHANNEL-AVERAGE `fbar` over smooth, not the four-sided `f` that
  `evaluate_rib` returns -- read the other way the model sits 1.8x off and
  no Re reaches the data. The smooth references are Han's lab's own, printed
  in NASA CR-4015 Eqs. (4)/(5) (Blasius, Dittus-Boelter); Han and Zhang's
  paper is not on disk, and swapping `f_s` for `0.046 Re^-0.2` moves the
  score 0.3%. Declared in `conventions.py`.

  **The geometry is inferred and the score moves with it.** The text states
  only Re 15,000-80,000; e/D 0.0625, P/e 10, square channel are borrowed
  from 91-GT-3 (same lab, same Re range). Recovered Re lands at
  34,000-65,000, inside the stated range -- a check that could have failed.
  Across Han's own e/D band 0.047-0.078 the bias runs -10% to +6%.

  The 60 deg V-shaped series stays unscored: no set expresses rib shape,
  and scoring it as a parallel rib would alias one shape onto another
  (#434). Both figure 4.53 series drop their `uncertainty: 0.08` to `null`:
  Han and Zhang's own band is unread, and a band borrowed from another
  paper is not a measurement (same treatment as fig4.193c).

- **`han_1989_narrow_channel`: Eq. 4.19 as its own selectable `RibCorrelationSet`,
  extending the Han rib family to narrow-aspect-ratio channels (`W/H` 1/4-1,
  `alpha` 30-90 deg).** Closes #402, the largest single block of unscored
  rib data left -- 37 previously-"not scored by any set" points in
  `fig4.48_{R,G}_WH{0.25,0.5}` now score.

  `G` is text-extracted: a hard ~20% switch in its leading constant at
  `alpha == 90 deg`, plus a further `W/H` correction below `W/H = 1/2` to
  both the constant and the `e+` exponent itself (`GShapeModel::
  NarrowChannelAlphaSwitch`). `R` has **no printed equation** below
  `W/H = 1` -- only a drawn figure line -- so its coefficients
  (`RAlphaShape::QuadraticAlphaTwoBand`, one quadratic-in-alpha per `W/H`
  sub-band) are this project's own fit to that line, and the set's
  `provenance` is `Fitted` rather than `Extracted` as a result. See
  `validation/cooling/extractions/han_ribbed.md`, item 40 and decision D6.

  Both new schema shapes are additive, following the same pattern #334
  used for `han_park_1988_angled`: `han_1988_orthogonal`,
  `rallabandi_2009_high_re` and `han_park_1988_angled` are provably
  unaffected, their full test suites passing unmodified. The `*_WH1_ref3`
  drawn-line and error-bar series stay deliberately unscored.

- **Taslim & Spring (1987), AIAA-87-2009 digitised -- the first rib source in
  this dataset that is not Texas A&M.** Northeastern plus General Electric.
  38 series over five figures: friction and Nusselt for eight two-side
  turbulated configurations and their one-side counterparts, 90 degree
  transverse in-line ribs, `P/e` = 10, Re 20k-190k, three aspect ratios.

  **It covers two things nothing else here does.** Taslim's AR 3.5 is Han's
  `W/H` = 0.286, inside the `1/4 <= W/H < 1/2` branch of Eq. 4.19 that #402
  will implement from a TAMU source. And `e/D` 0.083 at a square channel is
  the closest any independent measurement gets to Han's 0.047-0.078, which
  is what #403 item 2 is blocked on.

  **The aspect ratio is INVERTED relative to Han** -- Taslim's is
  height/width, Han's `W/H` is width/height -- and `e/D_H` is *coupled* to
  it, because `e` is fixed at three sizes and `D_H` follows. So low blockage
  occurs only at high aspect ratio, and no configuration here reaches Han's
  band at a square channel. A property of the rig: every independent
  arbitration of #403 item 2 will carry an extrapolation cost.

  **Five figures, each cross-checking another.** Calibration 0.6-1.8% on x,
  0.20-0.57% on y. The paper's own smooth-duct curve reproduces
  `0.023 Re^0.8 Pr^0.4` to **+0.30%** (fig4) and -1.77% (fig5), confirming
  the `Nu x 1e-3` ordinate without being told. Figures 4+5 reconstruct
  figure 9 to **5.4% MAE, -0.4% bias**. The paper's `Nu_T ~ Re^0.6` comes
  back out of the marks at 0.58-0.65.

  **Figure 12 misdraws one marker and it is the most valuable one.** Six of
  seven land within the x bow; AR 1.0 at `e/D` 0.083 is off by +0.047 in x
  and reads `f` = 0.0746 against figure 11's **0.0468** from eight marks --
  wrong in both coordinates, confirmed at 450 dpi. Figure 12's own fit curve
  gives 0.0437 there, so the curve is right and the symbol is misplaced.
  **Figure 11 supersedes figure 12 wherever they overlap.**

  **Figure 9's exponent is 0.2 and the figure appears to say 2** -- the
  superscript's decimal point does not survive the scan. The body text says
  ".2" twice; exponent 2 would put the ordinate at 35-523 on an axis ending
  at 5. Same class as the lost decimal in Han's `(e/D)^0.014`, worth 22%.

  **Committed UNSCORED.** The source publishes a passage-average *Fanning*
  friction factor and a turbulated-surface Nusselt number, not `R` and `G`.
  The conversion needs this family's own decomposition, whose `R` carries an
  aspect-ratio term inside the logarithm and is therefore a different
  quantity -- declaring that before applying it is exactly what #433 built.
  Runner work belongs with #402 and #403.

  Record: `validation/cooling/extractions/taslim_spring_1987_aspect_ratio.md`.

- **The scorecard answers #389's question 1 for the first time: a fidelity
  count per correlation set.** Fidelity was "partly built and entirely
  invisible" -- the checks existed, scattered across unit tests and
  extraction records, and nothing reported how many printed quantities an
  implementation had been checked against.

  ```
  Fidelity: checks against quantities the source PRINTS
    andrews_1986_effusion_internal             7
    baldauf_2002_sellers                       2
    han_1988_orthogonal                        1
    han_park_1988_angled                       1  (1 from figure cards)
    No printed-quantity check at all: florschuetz_1981_inline,
    mcgreehan_schotsch_1988_cd, mcgreehan_schotsch_1988_crossflow_cd,
    rallabandi_2009_high_re
  ```

  **A check compares against something the source PRINTS** -- a constant,
  a table entry, an equation the author evaluates himself, or an identity
  between two printed equations. Scoring against a digitised figure is
  not one: that is the basis the scorecard already reports, and it cannot
  separate "we transcribed it wrong" from "the model misses". Andrews
  Eq. (19)'s printed `0.27` against our `0.26983` can.

  **#389 expected this to be "mostly aggregation of checks that already
  run", and measuring showed it is not.** `verify.py` finds nine
  `printed-curve`/`printed-exponent` findings across the whole dataset
  and only ONE lands on a scored set -- the rest sit on `scores: null`
  correlation curves belonging to no set. So each substantial check is
  declared with its quantity, location, agreement and the test that pins
  it, and `test_every_declared_check_exists` resolves every named pytest
  node by collection and every gtest against the C++ sources.

  **The count is evidence, not a score.** `baldauf_2002_sellers` carries
  a check whose result is that Eq. (31) disagrees with the paper's own
  Table 4 by 36% -- fidelity evidence of the most useful kind, and it
  would be perverse for it to read as a defect. Sets with no check are
  NAMED, because "no check" and "nobody looked" both read as zero.

  The no-check list is derived, not declared: a hard-coded one reported
  `han_park_1988_angled` as having none while the machine count found
  one, contradicting itself in adjacent lines of the same report.

- **Measurement conventions are declared and checked, so scoring can no
  longer compare two different quantities** (#389's third remaining item).
  `validation/cooling/conventions.py`: a series declares what its y
  quantity is referenced to, a correlation set declares what it produces,
  and scoring matches them, applies a **registered** conversion, or
  refuses. No range check can see this failure -- every coordinate is in
  range, the quantity is simply not the same one.

  The project had hit it three times, and all three are now named in one
  place with their cost: Rohde's `Cd` on duct **total** against
  McGreehan-Schotsch's on **static** (1.33x at VHR 2); Andrews
  86-GT-225's `h_m` over the **hole length** against 88-GT-290's `h` on
  the **plate area** (-60% versus -10%, #430); and a gas-side `h` quoted
  without its unblown baseline (1.75x, #431).

  The whole dataset resolves: **100 direct, 3 converted, 0 undeclared**,
  and `test_every_scored_series_declares_a_convention` makes a new scored
  source fail until its convention is written down.

  **The declaration is a table keyed (source, quantity), not a field on
  103 series.** A convention is a property of how a rig was instrumented,
  not of the individual curve; one auditable table beats the same string
  copied 103 times and drifting. `SERIES_OVERRIDE` covers the exception.

  **`orifice_runner` now applies the Rohde conversion because the registry
  says to, not because an axis is named `velocity_head_ratio`.** The
  conversion was always correct; what was missing is that nothing would
  have noticed a series arriving already on the static basis. It now
  raises rather than converting twice -- 41% high at VHR 2.

  **A registered conversion must be geometry or algebra, never a fit.**
  Han's `G_bar = 1.2 G` is the case that tests the rule: same shape, and
  deliberately NOT registered, because the measured ratio spans
  1.096-1.413 with 69% of the variance between rigs (#401). It stays a
  published constant applied with its cost measured.

  **`eta_adiabatic` and `eta_overall` are declared incompatible**, with
  the reason recorded: an adiabatic wall passes no heat so its
  effectiveness carries no coefficient, and an overall-eta correlation
  used as a closure would double-count the convection the resistance
  network computes (#387).

  Falsification caught a promise with nothing behind it: the module
  docstring named `test_the_undeclared_count_does_not_grow` as the
  ratchet keeping `undeclared` from becoming permanent, and that test did
  not exist. Written, and it goes red when a source is dropped from the
  table.

- **The validation harness no longer pools in-domain with out-of-domain
  scoring, and a correlation's accuracy now carries where it came from**
  (#389, two of the four remaining items).

  **`AccuracyProvenance` and `StatedAccuracy`.** `#389` asked for the rib
  sets' `accuracy_R`/`accuracy_G` to be used as "the source's own band".
  Reading them showed that cannot be done as written: **only
  `han_1988_orthogonal`'s figures are an author claim.**
  `han_park_1988_angled`'s 10.5%/8.8% were measured by this project
  *through `evaluate_rib`* -- the very code path they would be judging --
  and `rallabandi_2009_high_re`'s 6.9% from the printed equation against
  digitised points. That is the same circularity #415 removed from the
  dataset's `uncertainty`, one level down and in C++.

  | set | accuracy_R | accuracy_G | usable as a band |
  |---|---|---|---|
  | `han_1988_orthogonal` | 6.0% Stated | 8.0% Stated | **yes** |
  | `han_park_1988_angled` | 10.5% Measured | 8.8% Measured | no |
  | `rallabandi_2009_high_re` | Unstated | 6.9% Measured | no |

  `usable_as_band()` is true only for `Stated`, and the report prints the
  provenance beside every figure so a reader cannot mistake one for the
  other. **`Unstated` now carries NaN, not 0.0** -- the old sentinel read
  as "perfect agreement", which is the one thing an absent measurement
  must never look like. `validate_rib_set` refuses a value without a
  provenance or a provenance without a value.

  **In-domain and out-of-domain are partitioned, never pooled.** Every
  point already carried `extrapolated`; the scorecard totalled it into a
  column beside an error that mixed both regimes. `build` and `rollup`
  now split on it and tag the out-of-domain rows `{out}`. What that
  exposes:

  - `baldauf_2002_sellers` scores 138 points of which **all 138** are
    outside its envelope, so its 48% MAE is an extrapolation figure and
    now says so.
  - `han_1988_orthogonal`'s largest fidelity row splits 46 in-domain at
    6.1% MAE against 52 out-of-domain at **4.6%** -- extrapolation
    scoring BETTER. Han's `e+ >= 50` floor and 60,000 Reynolds ceiling
    are conservative there, so filtering out-of-domain points would have
    discarded the better half. The partition is a report, not a filter,
    and that case is pinned by a test.
  - One of Lau's nine `G` points needs a Reynolds number above Han's
    stated ceiling, and both `fig4.47` panels straddle the box -- so the
    `measured` accuracy figures above, taken from those very panels, mix
    the two regimes. One more reason they are not bands.

  Falsification found a real gap and it is worth recording: the rollup
  test passed while the rollup stopped keying on domain at all, because
  it checked the label and the domain field against each other rather
  than against the data. It now asserts every rollup row is PURE --
  `n_extrapolated` is 0 or equal to `n` -- which is the invariant that
  binds the bucketing.

- **Andrews 88-GT-290 Figure 10 digitised, and #387's overall-effectiveness
  output shipped -- with a measured negative result attached.**
  `EffusionPlateElement.overall_effectiveness(state_in, h_gas, T_gas,
  U_gas=None)` returns `eta = (Tg - Tw)/(Tg - Tc)` and the wall
  temperature as OUTPUTS, with the gas side as the caller's input. Scored
  against Andrews' two effusion plates:

  | plate | D mm | X/D | n | bias | MAE |
  |---|---|---|---|---|---|
  | C | 3.27 | 4.7 | 72 | **+3.1%** | 4.6% |
  | B | 2.16 | 7.1 | 66 | **+23.7%** | 23.7% |

  Reported as two rollup rows, never pooled -- a pooled +13% would
  describe neither. Grouped per PLATE and not per jet regime: a
  velocity-ratio split looks informative (+7% below VR=1 against +19%
  above) but both plates cross VR=1 inside their own G range, so it
  reports the plate mix as physics. Holding the regime fixed, plate B is
  +21.2% where plate C is +2.2%.

  **The closure carries no separate film term, and that is not because
  there is no film.** An external film changes both the driving
  temperature and the coefficient -- `q = h_f (T_aw - T_w)` with
  `T_aw = Tg - eta_f (Tg - Tc)` -- and an adiabatic effectiveness
  correlation supplies only `eta_f`, because an adiabatic wall passes no
  heat and so measures no coefficient. Offered alone at the smooth-duct
  `h_g`, Andrews' data admits at most `eta_f = 0.105` (and NEGATIVE film
  at all 66 of plate B's points) against the 0.27-0.58 Baldauf gives over
  ten rows. Admitted with its matching augmentation, the same data demands
  1.6-4.3x the smooth-duct coefficient. Only the pair is identifiable from
  one eta curve, so `h_gas` is the pair -- a form with zero free
  parameters, where the split form fits anything.

  **What Baldauf does contain, corrected.** The counter-rotating vortex
  pair entraining hot gas under a lifted jet changes `T_aw`, which is
  exactly what IR on an insulated wall records -- so Baldauf HAS it, in
  Eq. (38), decaying as steeply as `mu^-4.57` at plate B's spacing.
  Andrews' "jet stirring of the cooling film" sentence describes that same
  film destruction, so it names a mechanism the correlation models rather
  than a missing one. An earlier draft of the provenance record had this
  backwards.

  **The finding: a gas-side residual no film correlation can close.**
  Plate B needs about 1.7x plate C's gas-side coefficient, and the
  smooth-duct value divides out of that ratio, so it rests on no
  assumption about the rig. Scaling Baldauf's `eta_f` from 0 to 1.25x
  moves the absolute requirement a long way (plate B's from 2.07 to 5.99)
  and the ratio only from 1.96 to 1.60 -- **so the split is bounded away
  from 1 whatever film model is chosen.** Closing it needs a gas-side
  augmentation `h_f/h_0` for full-coverage effusion, which is a heated-wall
  measurement; a better lift-off film correlation would not do it.

  No threshold is offered and none is applied: the requirement collapses
  on none of the velocity ratio, the blowing ratio or the momentum flux
  ratio. `jet_regime()` and the element report all three beside the
  result so a caller can see where they are.

  Six further Fig. 10 curves are committed unscored -- A/B, A/C and
  impingement A need a closure that does not exist here, and the porous
  wall, Lamilloy and Transply curves are other technologies. They earn
  their place: **the paper's own Table 3 tabulates the ratios BETWEEN
  these curves**, and reading it as `(eta_num - eta_den)/eta_den`
  reproduces all nine evaluable cells to 0.57 percentage points -- a third
  channel that validates the curve identification, not just the
  coordinates.

  **The page is rotated and the raw picks are committed anyway.** The
  frame corners give 0.94% skew one way and 0.96% the other -- one real
  rotation, measured twice. But the digitiser's affine already absorbed
  it: the four-corner inverse makes the y ticks *worse*, 0.00178 to
  0.00255. This is the Rohde Fig. 10 check run again and coming out the
  other way.

  **The gas side is TWO inputs, not one, because a film is two numbers.**
  `overall_effectiveness(state_in, h_gas_unblown, T_gas, U_gas=None,
  gas_augmentation=1.0, eta_film=0.0)`. The standard form is
  `q = h_f (T_aw - T_w)` with `T_aw = Tg - eta_f (Tg - Tc)`, so a film
  owes both a driving temperature and a conductance `h_f/h_0`; an
  adiabatic correlation supplies only the first. A single `h_gas`
  argument silently meant "the augmented one" and a caller would have
  passed a clean-wall value -- a factor-of-two error with no symptom.
  Both new arguments default to their no-film values, so omitting the
  pair is consistent rather than half-right.

  **`EffusionPlateElement` gains `internal_Nu_multiplier`**, matching
  `Nu_multiplier` on `ConvectiveSurface`. The element previously had no
  coolant-side tuner at all -- exposed by verifying that the knobs
  actually reach a target they are offered for.

  **Both knobs verified to reach, and one of them recorded as the wrong
  dial.** On plate B at G = 0.6 (measured 0.607, untuned 0.761, +25.3%):
  `gas_augmentation = 2.060` closes it exactly, as does the physical pair
  `eta_film = 0.318` with `gas_augmentation = 4.322`. So does
  `internal_Nu_multiplier = 0.486` -- by HALVING the coolant side, which
  is the wrong direction, since that correlation runs 10.4% *low* against
  Andrews' own Fig. 8. Reaching the target is not the same as being
  attributable, and using the internal knob to absorb a gas-side error
  would be reward hacking with a user-facing dial.

  **And one scalar is a rig match, not a model.** Fitting
  `gas_augmentation` at G = 1.0 takes plate B from 23.7% MAE to 5.8% --
  a real improvement -- but the residual runs -28.7% to +3.5%, because
  the required augmentation falls at low G. The average closing hides the
  tail.

  **Augmentation ratios need their baseline stated.** The required
  `gas_augmentation` for this rig is 3.1-4.3 against a fully developed
  Dittus-Boelter and 1.6-2.5 against the same duct thermal-entry
  corrected -- a factor of 1.75 from that choice alone, which is the same
  class of definitional error as the `h_m`-vs-`h` trap in #430.

  Two figure-vs-text conflicts recorded, neither resolved: the plot prints
  `Tg = 744K, Tc = 293K` where the text says 750 K and 295 K, and `Z =
  6.4mm` where the text says "An 8 mm impingement gap was used
  throughout". Neither touches the scored curves.

- **Andrews 86-GT-225 Figure 8 digitised: 41 points, and the effusion
  internal correlations now report a FIDELITY row as well as an accuracy
  one.** Four test plates against the paper that publishes the
  correlations, at **-10.4% bias, 11.1% MAE**. The scorecard reports it
  separately from the `andrews1988` accuracy row (-13.5%, 10 points),
  because a fidelity miss and an accuracy miss are different findings.

  | plate | D mm | X mm | L/D | technique | bias |
  |---|---|---|---|---|---|
  | a | 1.178 | 15.20 | 5.38 | Drilled | -6.2% |
  | b | 0.642 | 6.10 | 9.92 | **Spark Eroded** | -19.6% |
  | c | 0.897 | 6.10 | 7.05 | Drilled | -5.5% |
  | d | 1.298 | 6.10 | 4.85 | Drilled | -12.5% |

  **The fidelity miss is NOT a transcription error, and that is the
  finding.** A 10% gap on the author's own data would normally mean our
  reading of the equations is wrong. Figure 10 of the same paper rules
  that out: it plots Andrews' own evaluation of his own Eq. (19) and
  combaero reproduces it to 1.5%, at a Reynolds number that sits inside
  Fig. 8's range (18 points below, 23 at or above). So the two figures of
  one paper separate the failure modes outright -- **our transcription is
  right to 1.5%, and Andrews' correlation under-predicts Andrews' own
  measurements by about 10%.** Neither figure alone supports that.

  **Two papers, two definitions of "h", worth a factor of up to 9.8.**
  86-GT-225 defines `h_m` "over the hole length" (the hole internal area,
  `pi D L`); 88-GT-290 defines `h` "based on the surface area, A" with
  `A = X^2 - pi D^2/4` (the plate area). Same group, same figure number,
  different quantity. Scored on the wrong convention these points read
  **-60%**. The two are now distinct `y_axis` values -- `h_hole_length`
  and `h_internal` -- so the runner cannot confuse them silently.

  **Plate b is the worst fit and the only spark-eroded plate.** Recorded
  as an observation and nothing more: one plate is not evidence, no
  correction is applied, and excluding it gives -8.2% over 33 points --
  a smaller number for a smaller dataset, not a better model.

  The residual is not explained by Reynolds number, L/D or pitch; each is
  checked and reported in
  `validation/cooling/extractions/andrews_effusion_internal_h.md`.

  **The axes are log and were digitised on a linear calibration** -- the
  third figure in this repo caught by a calibration file failing to
  reproduce round numbers. `fig8_calibration_as_read.csv` is committed
  UNCONVERTED as the evidence: the x ticks at true 0.5 and 1.0 read 1.1299
  and 1.5689 where a linear pick of a log axis predicts 1.1218 and 1.5619.
  The inverse transform recovers the two ticks it does not use to 1.3% and
  1.5%, and the paper's own statement that G ran "from 0.1 to 1.7 kg/sm2"
  confirms the converted span of 0.089 to 1.79.

- **Andrews 86-GT-225 Figure 10 digitised: an IMPLEMENTATION check on the
  effusion internal correlations.** The figure plots four series on one
  axis -- the measurement, and the author's own evaluation of three
  correlations including his Eq. (19), which is what combaero implements.
  That makes it a different and stronger question than `andrews1988`
  Fig. 8's model-vs-data check, and the two together separate the failure
  modes #389 exists to keep apart: **did we transcribe the equations
  correctly** against **do the equations describe a rig**.

  | series | checks | agreement |
  |---|---|---|
  | triangles, Eqs. (12)-(14) | `mills_entry_length_factor` | **0.8% mean, 1.6% worst** |
  | squares, Eq. (19) | `effusion_internal_nusselt` | **1.5% mean, 2.5% worst** |

  That is digitisation precision, and it means the -13.5% against Fig. 8 is
  a **model limitation rather than our bug** -- a distinction the earlier
  PR could state only as an inference.

  **Plate A is not an outlier, it is the only pitch test in the figure.**
  Its measured point sits at 4.09 where the Eq. (19) curve is 2.73, which
  looks like a 50% miss. The paper explains it: plate a is a 10 x 10 array
  on a 152 mm plate, so its pitch is 15.2 mm against 6.08 mm for the 25 x
  25 plates -- and Eq. (18)'s approach term scales with `X/(pi L)`.
  Evaluated at plate a's own pitch, our Eq. (19) gives 4.21 against the
  measured 4.09, agreeing to **2.8%**. The other three points share one X
  and cannot distinguish the pitch dependence at all.

  The paper's own words: "the smaller number of holes has increased the
  heat transfer coefficient at a fixed mass flow rate".

  Recorded with the coverage split the perturbations revealed: Fig. 10
  spans `w = D/L` of 0.10 to 0.21, so it pins Eq. (14)'s LINEAR term and
  not its high-order ones -- changing `w^3`'s 58.6 to 56.6 leaves it green.
  Those are pinned instead by the branch-junction test at `w = 0.5`.
  Neither test alone constrains the polynomial; together they do.

- **`EffusionPlateElement.internal_heat_transfer()`** wires the Andrews
  86-GT-225 correlations into the element, so a solved panel reports its
  coolant-side coefficient in `diagnostics()` rather than only through the
  validation harness. Returns `Re_hole`, `Nu_approach`, `Nu_throat` and the
  coefficient **both ways** -- `h_hole_area` on the hole internal surface
  and `h_plate_area` on the approach area, Andrews' own definition -- plus
  the `area_ratio` between them, 3.46 for his plate C. Reporting both is
  the point: using one where the other is meant is a factor-of-three error,
  and the identity `h_hole * A_hole == h_plate * A_approach` is asserted.

  Cross-checked against the runner scored on Andrews Fig. 8: the element
  and the harness agree to **0.02%**, the residual being the element's
  rounding to a whole number of holes.

  The external film is deliberately excluded -- Andrews shows internal and
  film are not additive in effectiveness, so an overall-effectiveness
  correlation must never be the closure here. No flow returns no
  coefficient rather than a zero.

  A falsification that found **nothing** is recorded as its own test:
  substituting `wall_thickness` for the drilled `hole_length` in the
  approach term passed every existing test, because they all used
  90-degree holes where the two are equal. At 30 degrees it doubles the
  approach Nusselt number silently.

- **Effusion plate INTERNAL heat transfer, closing #387's internal side.**
  Andrews 86-GT-225: a coolant hole cools its wall in two places, and the
  paper's point is that the approach dominates at the Reynolds numbers
  effusion runs at.

  - `mills_entry_length_factor` -- Mills' `R_Nu` for a sharp-edged entry,
    Eqs. (13) and (14)
  - `effusion_approach_nusselt` -- Sparrow's hole approach, Eq. (18)
  - `effusion_throat_nusselt` -- Mills' short-hole throat, Eq. (12)
  - `effusion_internal_nusselt` -- the two summed, Eq. (19)

  `Re` is on the hole diameter and every `Nu` on the **hole internal area**,
  which is what makes them summable. A plate-area coefficient needs
  `* k/D * A_h/A`; that factor is 3.46 for Andrews' plate C, so it is not
  optional, and the conversion is left to the caller because it needs only
  geometry.

  **Scored against Andrews' own Fig. 8: bias -13.5%, RMSE 14.7% over 10
  points, labelled `accuracy`.** Cross-source, not fidelity: the
  correlations are the 1986 paper and the data the 1988 one, and the
  validation policy counts same-lab-different-study as cross-source. The
  first `SET_ORIGIN` entry used the bare author name, matched the
  `andrews1988` source and reported fidelity -- corrected to name
  `86-GT-225`.

  The miss is **one-sided**, so it is systematic rather than scatter,
  though it sits inside the source's own +/-16% spread and the two worst
  points are the two already flagged as off-trend in the source. Ruled out
  as causes: coolant temperature (unstated; worth 4 points across
  280-320 K), the `G` definition (per gross plate area, confirmed by
  `N = 4306` against `1/X^2 = 4305.8`), and hand-typed properties -- a
  first spike read -8% using values closer to 400 K air, where combaero's
  own properties give -13.5%.

  The scan is poor, so each equation was checked against something a
  transcription could not fake: the two `R_Nu` branches are independent
  polynomials that must meet at `L/D = 2` and agree to **1.3e-3**; `R_Nu`
  must decay to 1 for a long hole; and Eq. (19)'s printed `0.27` must equal
  `0.881 X/(pi L)`, which gives **0.26983**.

  `validation/cooling/effusion_internal_runner.py` wires it into the
  scorecard, so `andrews1988` is scored rather than carried unused. Record
  in `validation/cooling/extractions/andrews_effusion_internal_h.md`,
  including one thing this repo previously asserted and the paper does not
  say: Eq. (19) is a plain sum of Nusselt numbers, not a temperature chain.

- **Gao's alpha coefficients ARE published, and this repo said they were
  not.** Section 4.3.1 of Gao et al. (2025): "The empirical coefficients in
  Equation (5), a and b, were determined to be 12 and 0.9465,
  respectively." `film_superposition.md` recorded that the paper never
  prints them and that the search had been exhaustive -- it had looked for
  `a = 12` while the paper states the value in prose. Now recorded in code
  as `film_superposition::gao_a_case1` / `gao_b_case1`, with the header,
  both API references and the Python docstring corrected.

- **They still are not defaults, and the reason is arithmetic rather than a
  fit.** Since `a r/(a r + 1) >= 0`, that pair gives `alpha >= b = 0.9465`
  for **every** `r`. Murray's 5.75D staggered plate needs a per-row alpha of
  0.850, 0.864 and 0.686 at M = 0.19, 0.48 and 0.96 -- all below the floor,
  so no choice of `r` reaches them. The published pair can damp a row by at
  most 5.35%; this plate needs 15% to 31%. The likely cause is structural:
  Eq. (5) carries no streamwise-spacing term, and Murray's rows are 3.7x
  closer than Gao's calibration plate, which is the variable Murray proves
  dominant.

- **A superposition correction fitted on Murray, and the decision to keep
  exactly one.** #425 made `murray2018` the first family whose curves all
  fall one side of the model. For one geometry `r` is fixed per series, so
  the fit is a single-parameter scan:

  | M | Sellers | best alpha | MAE | best C | MAE |
  |---|---|---|---|---|---|
  | 0.19 | 33.1% | 0.850 | 15.8% | 2.15 | 13.7% |
  | 0.48 | 43.6% | 0.864 | 21.4% | 2.16 | 17.1% |
  | 0.96 | 112.6% | 0.686 | 23.7% | 4.29 | 16.2% |

  At M = 0.96 one knob takes the error from 113% to ~16-24%, close to the
  paper's own 15% experimental uncertainty.

  **The older "third category" coupling form `eta = eta1 + eta2 -
  C eta1 eta2` was evaluated and rejected**, so combaero continues to ship
  Gao's alpha alone. It fits 2-7 points better and shares the identity case
  (`C = 1` is exactly Sellers, verified to twelve digits) -- but it is
  **not bounded**: a plain sweep over C, eta and row count reaches
  **-1.7e5**, and it carries a hard ceiling at `eta = 1/C`. Gao's alpha
  never leaves [0, 1] over the same sweep, by construction. A network solve
  transits odd states during Newton iteration, so a correction that
  diverges there cannot be the one that ships, and a fit advantage inside
  the measurement band is not worth buying with divergence. alpha is also
  derived rather than asserted, and being a per-row vector it already
  subsumes what the coupling family wanted C for.

  Neither set of fitted values is shipped as a default -- they are tuned to
  Murray's plate and no one else's. Full record, including two earlier
  misreadings of the paper that were corrected, in
  `validation/cooling/extractions/gao_alpha_fit_on_murray.md`.

- **Murray & Ireland (2018) Figure 6, and the decomposition it makes
  possible.** Five digitised curves on a 5.75D-pitch effusion plate: three
  PSP measurements (M = 0.19, 0.48, 0.96) and two of the paper's OWN Sellers
  superpositions over a single-hole CFD field (M = 0.19, 0.96). The
  superposition curves carry `kind: correlation` and `scores: null` -- they
  are model output, never scored against, but they are the reference that
  separates the two error sources #420 could not tell apart.

  **The superposition is the dominant error, not the closure.** Their
  per-row input is CFD and ours is Baldauf, so with Sellers common to both:

  | M | superposition (theirs/exp) | closure (ours/theirs) | total |
  |---|---|---|---|
  | 0.19 | 1.18x | 1.13x | 1.33x |
  | 0.96 | **1.74x** | 1.28x | 2.24x |

  The superposition error roughly quadruples its excess with blowing while
  the closure barely moves -- which is what justifies spending a correction
  on the superposition rather than on the per-row term.

  **And this family CAN carry Gao's alpha, unlike `andrei2014`.** Scored
  against Baldauf + Sellers the bias is +23.5%, +22.3% and +96.7% -- all
  over-predicted, because M 0.19-0.96 sits entirely below the lift-off peak
  that drove Andrei's BR 2-3 curves the other way. alpha only ever reduces,
  so a family that falls on one side is a prerequisite for fitting it.

  The digitisation reproduces the paper's own prose independently:
  superposition/experiment is 1.74x mean at M = 0.96 rising to 2.09x by
  x/D 28, against the paper's "around twice at a blowing ratio of
  approximately one", and only 1.18x at M = 0.19.

  Notable extraction findings, all recorded in the metadata header: the row
  comb is **10 rows at 2.875D, not 5 at 5.75D** -- the plate is staggered so
  the film meets a row every half pitch, measured from the curves' own
  injection points and confirmed by the hole count (5 primary + 5 staggered
  against the paper's 4 x 5 array of 40 holes); the source file named
  "M=0.4" is the **M = 0.48** panel; and calibration ticks digitised in data
  coordinates land within 0.32% (x) and 0.41% (y) of the printed round
  numbers across all three panels independently.

- **`b_0_override` on the Baldauf film correlation**, so the paper's Table 4
  reading of Eq. (31) can be tested. Both entry points take it after
  `status`; NaN (the default) uses Eq. (31) as printed, and Python takes
  `b_0=None`.

  **Deliberately an override, not a second formula.** The paper contains no
  alternative: Table 4 states one number at the Table 3 conditions, and
  reaching it needs `sin(...) = -0.167` where the printed equation gives
  `+0.470` -- no sign, unit or angle convention makes that argument negative
  at 30 deg. Scaling Eq. (31) by the 0.737 ratio would be an invented
  correction, so the choice stays with whoever has evidence for it. Also
  confirmed dead: the "treat the sine argument as degrees" reading is
  *algebraically identical* to the printed one, because the 0.86 scale
  passes straight through a linear argument.

- **Corrected: Eq. (31)'s discrepancy is governed by `x/D`, not by `M`.**
  The header claimed "under 5% below M ~ 0.5 and up to 50% at M = 2.5". That
  misattributes it -- `b_0` reaches the model only through
  `b_1 = b_0/(1 + M^-3)` (Eq. 32), the *descending-branch* gradient, so near
  the hole it does nothing whatever the blowing. Measured with the new
  override across the whole envelope:

  | | effect |
  |---|---|
  | `x/D <= 20`, any `M` | under 1% |
  | `x/D = 200`, `M = 2.0` | -19% |
  | `x/D = 400`, `M = 2.0` | -28% |
  | worst in the box (`s/D` 5, 90 deg, `Tu` 0.0035, `M` 2.5, `x/D` 400) | **-75%** |

  Practical consequence: at effusion row spacings -- Murray's 5.75 D,
  Andrei's 9.15 D -- the two readings differ by **under 1%**, so no effusion
  dataset can arbitrate Eq. (31), and it cannot explain the Baldauf +
  Sellers failure measured in #420. Only far-downstream single-row data
  could settle it. Pinned by
  `EquationThirtyOneMattersFarDownstreamNotAtHighM`.

- **Multi-row film/effusion scoring, and the first verdict on the Baldauf +
  Sellers chain.** `validation/cooling/film_runner.py` scores superposed
  film effectiveness against adiabatic effusion-plate data, wired into the
  scorecard as `baldauf_2002_sellers`. Two datasets land with it:
  `andrei2014` (six adiabatic curves, BR 1-3 x DR 1.0/1.5, extracted exactly
  from the publisher's PDF vector stream) and `andrews1988` (effusion
  internal heat transfer, Fig. 8 plate C).

  Both sides are reduced to per-row means over each pitch interval, and the
  prediction is averaged over **the same abscissae** as the measurement --
  comparing an interval mean against a point value at the interval centre
  overstates the first row badly, where eta climbs from zero.

  **The chain does not describe this plate, and the runner says so.** Pooled
  over 108 row-means: MAE 45.5%, RMSE 53.1%, every point flagged
  extrapolated, basis reported as `accuracy` (Baldauf's own data is not in
  this repo, so `fidelity` is unreachable by construction). The sign of the
  error flips with blowing rate:

  | | DR 1.0 | DR 1.5 |
  |---|---|---|
  | BR 1 | +61.3% | +31.9% |
  | BR 2 | -7.8% | +17.0% |
  | BR 3 | -80.9% | -68.4% |

  The cause is jet lift-off. Baldauf's single-row effectiveness peaks near
  M ~ 0.7-1.0 at any spacing of s/D 3 or wider and collapses above it -- at
  Andrei's s/D it is 0.081 at M = 1 and 0.007 at M = 3. A dense 18-row plate
  does not behave that way; the measured effectiveness keeps *rising* with
  blowing. This is a single-row correlation failing to transfer to a compact
  array, not an implementation error, and it is not an artefact of
  extrapolating s/D: the early peak is already there at s/D = 3, inside the
  published envelope.

  **Gao's per-row alpha cannot be fitted on this family**, and a test now
  records why. alpha is bounded [0, 1] and only ever reduces the superposed
  value, but three of the six curves are over-predicted and three
  under-predicted, so no single (a, b) can serve both. The test goes red if
  a future change puts every curve on one side -- which is exactly when the
  original plan becomes viable again.

- **Multi-row film superposition** -- `film_superposition_sellers`
  (Gao Eq. 1), `film_superposition_corrected` (Gao Eq. 7, Sellers with a
  per-row mainstream temperature correction), `mainstream_temperature_correction`
  (Eq. 5's form), `equivalent_slot_width` (Eq. 9) and
  `equivalent_blowing_ratio` (Eq. 10), each with `..._and_gradient`
  variants where a solver needs them.

  **Sellers overestimates, and worsens as rows accumulate** -- Gao et al.
  (2025) measure exactly that. It is worst where effusion lives, many
  closely spaced rows, which is why #386 and #387 were designed together
  rather than in sequence.

  The correction `alpha` is derived, not invented: Gao's Eqs. (3)-(4) are an
  energy balance on mainstream entrained into the boundary layer at each
  injection, so `alpha` is the fraction of the film's temperature deficit
  surviving to the next row. It is bounded in [0, 1], and `alpha = 1`
  reproduces Sellers exactly. This replaced an earlier proposal to tune with
  an invented coverage exponent.

  **Gao publishes alpha's form but not its fitted coefficients**, so `a` and
  `b` are required arguments rather than defaulted -- a made-up default
  would acquire an authority it has not earned. Fitting them against Andrei
  (2014) and Murray (2018), labelled tuned, is the next step.

  Gradients are built from explicit partial products rather than by dividing
  the total, so a fully effective row does not produce a NaN.

- **Baldauf et al. (2002) film-cooling effectiveness** --
  `film_effectiveness_baldauf_2002`, the laterally averaged adiabatic
  effectiveness downstream of one row of cylindrical, streamwise-inclined
  holes. Forty equations, plus `..._and_derivatives` giving
  `(eta, deta/dM, deta/dP)` analytically through forward-mode dual numbers
  over the same templated equation chain, so value and derivative cannot
  drift apart.

  Chosen because the alternatives fail where effusion and film cooling
  actually live: the review by Xia et al. (2024) notes that single-hole
  models are "completely wrong in the near wake region" and older row
  correlations "give an unrealistic maximum at the ejection position". This
  one is valid from the ejection point and carries the **adjacent jet
  interaction** -- the lateral hole-spacing effect driving jet lift-off -- as
  a correlated parameter. `deta/dM` correctly changes sign through the
  lift-off peak.

  `eta` is `(T_G - T_AW)/(T_G - T_C)`, already combaero's convention, so it
  composes with `adiabatic_wall_temperature` unchanged.

  **Angle units were measured, not assumed.** The API takes degrees; the
  paper's trigonometry is in radians. Seven Table 4 coefficients spanning six
  trig forms reproduce on radians and fail on degrees by 10-30% -- though one
  of them, `b*_T`, differs by only 0.06% and would have passed a units bug
  silently.

  **Eq. (31) is implemented as printed and contradicts the paper's own
  Table 4** by 36% in `b_0`, while the other 17 checkable coefficients
  reproduce to 2e-6. The printed form was confirmed from four independent
  channels; no angle convention, unit choice or single-token variant reaches
  the table's value. Under 5% effect on `eta` below `M ~ 0.5`, up to 50% at
  `M = 2.5`. Both numbers are pinned in a test so the contradiction cannot be
  quietly rediscovered or quietly "fixed".

  Validated end to end against the paper's Fig. 14 (its own
  measurement-versus-correlation case) and Fig. 8's peak-effectiveness bands.
  Provenance in
  `validation/cooling/extractions/baldauf_2002_film_effectiveness.md`.

- **Lichtarowicz, Duggins and Markland (1965) long-orifice discharge
  coefficient** -- `DischargeCdCorrelation::Lichtarowicz1965`, Eqs. (7) and
  (12), replacing the member that had declared-and-refused since #409.
  Valid `l/d` 2-10 and **Re 10 to 2e4**.

  **This is the regime a cooling hole actually runs in.** Andrews' effusion
  plate C spans Re 432 to 8718 across its own measured data, entirely below
  McGreehan-Schotsch's `re_min = 1e4` floor, where that correlation holds Cd
  constant and reads **+21.7% high at Re = 432**. Cd genuinely varies 20%
  over that range.

  Refuses below `l/d = 1.5` rather than clamping: the source's own design
  recommendation is to avoid that geometry because Cd "varies rapidly" there
  and there is "the possibility of hysteresis", which a single-valued
  correlation cannot represent. `l/d` is held at 10 above the range, because
  Eq. (7) is linear and extrapolates to zero near `l/d = 97`.

  Closed-form, so `dCd/dRe` is differentiated analytically rather than
  through dual numbers; matches central differences to ~1e-9 relative.
  Provenance in
  `validation/cooling/extractions/lichtarowicz_1965_long_orifice.md`.


- **Validation policy: fidelity, accuracy and tuning are three questions,
  never one number** ([docs/VALIDATION_POLICY.md](docs/VALIDATION_POLICY.md),
  and a Hard Rule in `CLAUDE.md` -- it applies to the whole `validation/`
  tree, not only cooling).

  **Fidelity** -- does the implementation mirror the paper -- is judged
  **only on the author's own data**; a miss is our transcription.
  **Accuracy** -- what the model does on another rig -- is judged
  **cross-source and labelled cross-source in the report**; a miss is the
  model's limitation. Same lab, different study counts as cross-source.
  Matching a particular rig stays the user's job with
  `Nu_multiplier`/`f_multiplier`, and the harness never scores,
  recommends or applies a tuner.

  The scorecard renders `basis` per row and the rollup segregates on it, so
  the two can no longer be averaged together.

  **What separating them showed.** `han_1988_orthogonal`'s fidelity rests on
  **eight cleanly-sampled points** -- its own data survives only as an
  overplotted scatter -- while its cross-source accuracy reads *better*
  (78.3% within against 62.5%) purely because the cross-source data is
  tabulated. Pooled, that read as a single unremarkable number.

- **Absolute `R` is scoreable, unblocking six committed series** (#339
  plan 1a). `runner.py`'s `e+` path resolved `y_axis` to `G`, `G_bar` or
  `R_normalised` only, so sources that TABULATE the roughness functions --
  which print `R` itself, not the normalised ordinate figure 4.46 plots --
  had to be committed `scores: null`. The measurements were never in doubt;
  only the branch was missing.

  **Reviewer item L2 becomes a number rather than an argument.** Lau's `R`
  now scores MAE 11.3%, bias -11.3%, with **0 of 9 points inside Lau's own
  stated +/-10.9% friction band**. Two labs ~12% apart in ribbed-wall
  friction for nominally the same 90 deg configuration, reported on every
  run instead of argued in a document.

  **Eq. 4.17 is confirmed out-of-sample**: NASA CR-3837's `R` scores 100%
  within band at 90 deg and 83-86% at 75/60/45, inside the source's own
  6.6% friction uncertainty. Only 30 deg is poor (57%), which is also where
  Eq. 4.18 reads furthest low.

- **NASA CR-3837 (Han, Park and Lei, 1984) as a tabulated primary source**
  (#391). The first source in the cooling dataset whose every run is
  tabulated rather than plotted, and the only one printing the ribbed- and
  smooth-wall Nusselt numbers separately, so the ribbed-wall Stanton number
  is read rather than assumed. 33 runs across five rib angles.

  Scored **out-of-sample** against `han_park_1988_angled`: Eq. 4.17/4.18
  were fitted to `e/D` = 0.047 and 0.078 on a different contract, and this
  is `e/D` = 0.063, a blockage they never measured. `R` is confirmed (bias
  -2.2%, MAE 6.8%, inside the source's own 6.6% friction uncertainty); `G`
  reads low at every angle, which takes the figure's side in item 23 of
  `han_ribbed.md`.

  Deliberately **not** scored against `han_1988_orthogonal`: figure 4.46's
  right-hand box prints "Han (1984), 8,000 <= Re <= 80,000" at `e/D` 0.063
  for both `P/e`, a four-way match to CR-3837's own Table 2, so those marks
  are these runs and that correlation was fitted to them.

  Values are derived to Han's four-sided friction basis and ribbed-wall
  Stanton number; the printed `Rbar` is 5.35 where the same measurement on
  Han's basis is 3.195, so committing the printed columns would have been a
  67% definitional error. Page 128 is two-channel verified.

- **The cooling scorecard segregates by sampling completeness, and recovers
  the marks a digitisation could not pick** (#393). A paper prints graphs,
  not tables, so where a figure overplots symbol classes only the spatially
  isolated marks can be digitised -- and those are the ones furthest from
  the cluster centre. Rows are no longer pooled across that divide.

  **What it revealed.** Pooled, `han_1988_orthogonal` reported 73.0% of
  points within the source's stated band, against Han's own printed claim
  of 95% -- a 22-point gap that read as a deficiency in the correlation.
  Segregated, its completely-sampled series reach **91.1% within, MAE 2.9%**
  and its partially-sampled ones **65.3%, MAE 5.3%**. The gap was the pick,
  not the model.

  `partial` rows are labelled in the rendered scorecard as an **upper bound**
  on the error rather than an estimate of it.

  **Recovery** (`validation/cooling/recovery.py`). A mark that could not be
  picked is not lost: where a figure stacks two panels over one abscissa, a
  mark visible in either proves the run exists in both, and its ordinate is
  bounded by the envelope of the marks that are visible nearby. That turns a
  missing mark into an interval observation -- 26 of them, envelopes 1.3% to
  13.9% wide.

  The envelope is drawn from other digitised marks, never from the
  correlation plus its stated band. The latter is tighter and tempting and
  would bound the observation with the model being scored. Recoveries whose
  envelope exceeds 15% are discarded rather than counted as weak evidence,
  which is why figure 4.46 yields 25 and figure 4.51 yields 1 -- 4.46's
  classes hug one correlation, 4.51's genuinely separate.

  Containment is reported as its own column (`held`, currently 48%) and
  deliberately **not** folded into MAE/RMSE: a point gives a signed error, an
  interval gives containment, and averaging the two would let a loose bound
  flatter a set that a tight point disagrees with.

- **Lau's `Gbar` committed as an independent check of Han's `G_bar = 1.2 G`**
  (#392). The relationship previously rested on a single printed label on
  figure 4.47. Lau's Eq. (8) measures the same quantity -- `St_avg` over two
  ribbed and two smooth walls -- so it is a like-for-like comparison by
  another lab: N=9, MAE 3.5%, RMSE 3.9%, bias -3.5%, 88.9% inside Lau's
  stated +/-5.8% Stanton uncertainty. Han's factor sits a few percent below
  Lau's, recorded rather than reconciled.

- **Lau, Kukreja and McMillin (1990) as validation data** (#385), the first rib
  source outside Han and a primary paper rather than a textbook reprint. Its
  Table 2 prints the roughness functions as closed forms, so the 90 deg
  baseline needed no digitising at all (`extraction: tabulated`,
  `confidence: exact`).

  **`han_1988_orthogonal`'s heat-transfer roughness function is independently
  confirmed**: Han prints `G = 3.7 (e+)^0.28`, Lau `4.218 (e+)^0.257` --
  different coefficient and exponent, agreeing to RMSE 1.4%, bias -0.1% over
  `e+ = 106-636`, with every point inside Lau's own stated +/-5.8% Stanton
  uncertainty. A different lab, rig and decade.

  **The friction roughness function does not agree**: Lau's `R` runs 12.4-13.1%
  above combaero's constant 3.2, which converts to two labs ~12% apart in
  ribbed-wall friction factor for nominally the same configuration.
  Definitional causes were checked and ruled out. Recorded, not reconciled;
  committed with `scores: null` because `runner.py`'s `e+` path has no
  absolute-`R` branch.

  **Reflection folding is falsified.** `han_ribbed.md` recorded
  `alpha -> 180 - alpha` as the bounded alternative if Eq. 4.17 were ever
  guarded past 90 deg. Lau measured both sides: V-120 differs from V-60 by
  15-23%, and the paper's own conclusion says reversal lowers heat transfer.
  The open item is closed -- no guard should be written on that basis.

  **Lau's `Gbar` is deliberately NOT committed.** It is a four-wall average;
  Han's `G_bar` is the Prandtl-normalised `G Pr^-0.57`, which
  `han_ribbed.md` item 10 explicitly resolved as "not a four-wall average".
  Their ratios to `G` differ by only 2-3%, so scoring one against the other
  would read as a confirmation -- and `runner.py` multiplies anything tagged
  `y_axis: G_bar` by Han's Prandtl factor automatically. A test asserts the
  trap stays shut.


- **Adiabatic expansion factor for the McGreehan-Schotsch orifice** (#382):
  `mcgreehan_schotsch_1988_expansion_factor` / `_expansion_orifice` /
  `_expansion_nozzle`, and `orifice::mcgreehan_schotsch::expansion_*` in C++.
  Eqs. (4)-(7) of the source, which #381 left unimplemented -- the paper's
  mass flow is `W = Cd * Y * A * sqrt(...)` and only `Cd` had been shipped.

  **`Y` is for the incompressible formulation only.** combaero's
  `regime='compressible'` already solves the isentropic nozzle exactly via
  `nozzle_flow`, and Eq. (5) reproduces that solve to **0.008%** over
  `0.6 <= S <= 0.99` -- so applying `Y` there would correct for
  compressibility twice (worth >8% at `S = 0.7`). That agreement, between a
  1988 closed form and combaero's own independent numerical solve with
  real-gas properties, is also the strongest check on the transcription.
  `orifice.h`'s module comment is corrected to say which route owns the
  expansion term.

  **Eq. (7)'s blend weight saturates smoothly rather than with a hard clamp.**
  The paper prints a bare ramp with no upper bound, and our own chain reaches
  `X = 1.45`, extrapolating past a nozzle. A hard clamp would put `dY/dCd` at
  *exactly zero* outside `Cd` in `[0.82, 0.94]` -- a reachable design point --
  plus a discontinuous jump of 0.734 at both knees, which would break the
  finite-difference agreement required of the `(f, J)` work in #383. The
  smoothed form costs 0.47% in `Y` at worst and its `eps = 0` limit is exactly
  the paper's clamp, so the exact form remains available. Continuity is tested
  by how the derivative change scales with sampling interval, not against a
  threshold: the hard clamp holds 0.73353 at every interval, the smoothed one
  halves as the interval halves.

  `Y` also saturates at the critical pressure ratio, below which the
  isentropic form predicts decreasing flow; the proxy mass flux was confirmed
  to peak exactly at `S*`. Decisions D10-D13 in
  `validation/cooling/extractions/orifice_discharge_coefficient.md`.

- **`validation/solver_smoothness.py`**: a reusable scan for the three ways a
  blend or clamp stalls a Newton solve -- FLOOR (derivative exactly zero over
  a run), KINK (derivative discontinuous), DIVERGENCE (derivative unbounded).
  They are told apart by how the derivative behaves as the sampling interval
  shrinks, which also separates all three from an honest steep gradient.
  `assert_smooth()` for tests, with `allow=` so a known, documented limitation
  narrows the check instead of switching it off; `python -m
  validation.solver_smoothness` for a report.

  Written because reading the code does not find these: the expansion factor
  above had two bounds smoothed deliberately and still shipped a hard
  `min(S, 1)` putting a kink at the no-flow boundary, which the scan caught.
  The detector is itself tested against hazards with known analytic answers in
  both directions -- it must fire on `max(x,0)`, a hard clamp, and `sqrt|x|`,
  and must stay quiet on `tanh(40x)` (steep is not broken), on `x**2`'s
  isolated stationary point, and on smooth saturation.

- **Jet-plate / cooling-hole discharge coefficient**: `mcgreehan_schotsch_1988_cd`
  (Python) and `orifice::Cd_McGreehanSchotsch` plus the per-equation stages in
  `orifice::mcgreehan_schotsch` (C++), implementing McGreehan and Schotsch
  (1988), ASME J. Turbomachinery 110(2), 213-217, Eqs. (8)-(17). Closes the
  gap tracked by #375: a hole discharging plenum-to-plenum over Reynolds
  number, inlet corner radius `r/d`, orifice length `L/d` and inlet crossflow
  `U1/Vi`. The existing `Cd_sharp_thin_plate` family stays what it was, ISO
  5167 metering for a hole in a pipe run, which does not apply here.

  Validated against the source's own figures and against a second paper: the
  chain reproduces the paper's separately stated 0.60 baseline at
  `Re = 3.2e4` to 0.02%, lands on its nozzle asymptote to 0.18% when the
  corner-radius fit is extrapolated to the `r/d = 0.82` knee, and tracks the
  drawn curves of Figs. 3, 4 and 7 to 0.1-3%. For a bare plenum-fed jet plate
  (`r/d = 0`, `t/d = 1`, `Re = 1e4`) it gives 0.785 against Florschuetz,
  Truman and Metzger's independently recommended 0.79, and stays inside their
  measured 0.73-0.85 band across `t/d = 1-3`.

  Two behaviours are faithful to the source and deliberately kept: `Cd` is
  **not monotonic** in crossflow (it rises ~5.7% near `U1/Vi ~ 0.09` before
  falling, as drawn in Fig. 4 and supported by Rohde's NASA TN D-5467), and
  it is **held constant below `Re = 1e4`**, the stated validity floor, because
  Eq. (8) diverges below it. `U1/Vi` is a SUPPLY-side ratio and must not be
  fed Florschuetz's discharge-side `Gc/Gj`.

  `mcgreehan_schotsch_1988_crossflow_cd` applies Eq. (17) alone to a baseline
  the caller supplies -- for a plate whose zero-crossflow `Cd` is measured or
  taken from literature (Florschuetz's Table 1, say). This is also how the
  source uses Eq. (17) in its own validation.

  **Scored against Rohde's data** (NASA TN D-5467), as replotted in the
  source's Fig. 6 with the conversion to static parameters done by the paper's
  own authors: `t/d = 0.51`, 11 digitised points, bias `+5.95%`, RMS `7.54%`.
  The correlation reads high against Rohde and increasingly so with crossflow,
  which is the paper's own observation that Rohde's "basic values are lower".
  Quantified here for the first time: the chain's baseline exceeds Rohde's
  measured basic value by `+4.6%` at `(r/d, L/d) = (0, 0.51)`, `+8.6%` at
  `(0, 4.0)` and `+12.5%` at `(0.49, 1.06)`. Digitised data under
  `validation/cooling/data/mcgreehan_schotsch1988/`, scored by
  `validation/cooling/orifice_runner.py`. A 1.25x scaling of `U1/Vi` would cut
  RMS to 1.88% and is deliberately NOT applied: Eq. (17) was fitted to Rohde,
  Grimm, and Meyfarth and Shine pooled, so refitting it to the single series it
  is judged against would be tuning against the score.

  **Also scored against Rohde Fig. 10** (`validation/cooling/data/rohde1969/`,
  30 points, three inlet-edge conditions at `t/d = 1.06`), the only data in the
  project that tests the corner-radius term. Two independent marker counts
  agreed exactly (8/10/12); x-tick calibration 0.59%, y-tick 0.00185 in `Cd`;
  the `circle > square > triangle` ordering invariant holds everywhere.

  This surfaced a **measured limit of the correlation**: it reads high, and
  increasingly so with supply-side crossflow, but NOT equally across inlet
  radii -- `+40.6%` for a sharp edge over `VHR 1-10` against `+15.4%` for
  `r/d = 0.488`. The total-to-static conversion is identical for all three
  series at a given velocity head ratio, so it cannot produce an error that is
  monotone in `r/d`: **rounding protects against crossflow degradation more
  than Eqs. (12)+(17) predict**. A conversion-free ratio test puts it at bias
  `-13.1%`, RMS `15.7%`, while agreeing to `+3.3%` at low crossflow -- so
  Eq. (12) itself is sound and the interaction is what is missing. At
  `VHR >= 25` all three series converge on `+4%` to `+9%`, matching the
  independent Fig. 6 result. Nothing was tuned; pinned by
  `test_rohde_reveals_the_crossflow_limit_of_the_rd_term`.

  Practical bound: a plenum-fed jet plate has `U1/Vi = 0` and is unaffected. A
  duct-fed orifice is good to roughly `+8%` while `U1/Vi <~ 0.2`, and degrades
  badly beyond `U1/Vi ~ 0.35`, worst for sharp edges.

  Provenance, the cross-checks, and two errata found in the paper are recorded
  in `validation/cooling/extractions/orifice_discharge_coefficient.md`.

- **`make_constant_correlation(double Cd)`** (C++, `orifice.h`): pins the
  discharge coefficient to an explicit value. `CdCorrelation::Constant` could
  previously be selected but its value could not be set --
  `make_correlation()` always returned the 0.61 default.

- **Jet impingement wired into the network solver and GUI**: `ImpingementModel`
  and `SingleJetImpingementModel` (`python/combaero/network/components.py`),
  completing #337's re-add. `ImpingementModel` models ONE spanwise row of a
  jet array -- there is no single "channel Nu" for a whole array the way
  there is for a smooth or ribbed duct, since downstream rows see
  progressively more crossflow than upstream ones. A full array is a chain
  of these elements, one per row.

  Mass flow is recovered from `ConvectiveSurface.area`, the same convention
  every other channel model already uses -- no new one introduced. The
  element treats `area` as the row's own footprint, divides by the hole
  spacing product to get its hole count, and derives this row's own jet
  Reynolds number from the upstream state's total mass flow through that
  area. Converting an array's total/mean flow into a genuinely row-dependent
  distribution (Florschuetz's own Eq. 7) is deliberately left to whoever
  assembles the chain, not decided inside the element.

  Pressure drop is deliberately NOT modelled by either element: the jet
  plate's own orifice loss belongs on a proper `OrificeElement` upstream
  (with a real discharge-coefficient correlation), not duplicated here where
  it risks double-counting. `f`/`dP` reported by these elements are the
  crossflow/target channel's own plain smooth-duct values (Han's book gives
  no impingement-specific friction correlation, unlike ribs' own `R`/`f`
  relationship), borrowed to compute a `T_aw` in place of Eq. 4.2/4.3's own
  recovery-factor model, which is not wired in either.

  Both `impingement_correlation.h` functions (`single_jet_impingement_nu`,
  `jet_array_impingement_nu`) gained an analytic-derivative sibling
  (`single_jet_impingement`, and a new `dNu_dRe_j` field) for the wall-
  coupling Jacobian, per the project's (f, J) rule -- both FD-verified in
  `tests/test_impingement_correlation.cpp`, and the mdot-to-h composition
  chain independently FD-verified again at the element level in
  `python/tests/test_impingement_channel_element.py`, since that chain is
  new code the correlation-level check does not cover.

  Wired into the GUI the same way ribbed was: `ImpingementModelData` and
  `SingleJetImpingementModelData` (`gui/backend/schemas.py`), mapped in
  `graph_builder.py`, and a reachable "Impingement (Jet Array)" /
  "Impingement (Single Jet)" choice in `SurfaceEnhancementInspector.tsx` --
  replacing the old, unreachable, pre-0.7.0 impingement UI block that used
  the crossflow-less field names.

  Surfaced (not caused) a pre-existing, already-documented gap while
  extending `tests/test_channel_element_jacobian_fd.py` to cover the new
  elements: `ChannelElement._ribbed_residuals`'s own docstring already
  states its pressure-drop Jacobian has no temperature column, by design --
  adding "ribbed" to that test's surface matrix (which the file's own
  docstring had earmarked for #334, but was never done) required a
  documented skip for that one known combination, not a code change.

- **Jet impingement correlations re-added, with provenance** (`include/impingement_correlation.h`),
  closing issue #337. Two independent regimes:

  `goldstein_1986_single_jet()` -- a single free round jet on a flat plate
  (Goldstein, Behbahani and Heppelmann 1986, Han's Eq. 4.1). Reproduces the
  source's own closed-form check point (`Re=25000, R/D=5, L/D=7.75` ->
  `Nu=60` constant-heat-flux, `56` constant-wall-temperature) to within
  rounding.

  `florschuetz_1981_inline()`/`florschuetz_1981_staggered()` -- a jet array
  fed from a common plenum, where every row's Nu is degraded by crossflow
  accumulated from every row upstream of it (Florschuetz, Truman and Metzger
  1981, Han's Eq. 4.9/Table 4.1). The crossflow-to-jet mass flux ratio
  `Gc/Gj`, the correlation's defining input, is derived in closed form
  (`crossflow_to_jet_ratio_at_row`/`_at_x`, the paper's own Eq. 8) rather
  than left for a caller to supply -- it depends on `(yn/d)(z/d)` only, not
  `xn/d`, exactly as the source states.

  The prior removed implementation (#332) cited this same 1981 paper but
  took no crossflow-ratio input at all, making it structurally unable to be
  Florschuetz's correlation. This re-add was obtained and checked directly
  against a scan of the primary paper (not just Han's secondary reprint),
  page-image by page-image since it is a scanned copy -- every equation and
  Table 4.1's coefficients matched digit-for-digit, and the paper closed
  three gaps Han's reprint left open: the real validity range (Han's Table
  4.1 carries a misprinted rib-validity box instead, unrelated to jet
  arrays), the `Gc/Gj` closed form itself, and Eq. 4.7/4.8's attribution
  (the same paper's own less-tight alternate form, not a separate later
  one). See `validation/cooling/extractions/han_impingement.md`.

  Not yet wired into a network element or the GUI: converting an array's
  total/mean mass flow into a per-row local `Re_j` is a design question
  deferred to a follow-up, tracked under #337. Leading-edge/curved-surface
  impingement and the simpler forms (Eq. 4.6, Eq. 4.7/4.8) remain
  documented-but-unimplemented, matching how ribs staged their own re-add.

  Surfaced a related gap: combaero has no discharge-coefficient correlation
  for a jet-plate array (`orifice.h`'s `Cd_sharp_thin_plate` family is an
  ISO 5167 pipe-metering model and does not apply). Tracked separately,
  not blocking, as issue #375; the paper's own default (`C_D=0.79`) is used
  as-is.

- **Figure 6 digitised, validating `florschuetz_1981_inline()` against its
  own source**: all nine panels (three streamwise spacings times three
  spanwise spacings, three `z/d` symbol series each) -- 27 series, 242
  points, `validation/cooling/data/florschuetz1981/`.

  Scoring needs no Reynolds-number bisection, unlike every rib figure: the
  paper's own algebra makes `Nu/Nu1 = 1 - B[(z/d)(Gc/Gj)]^n` exactly, since
  `Nu1` is the same correlation's value at `Gc/Gj=0` and the
  `A*Re_j^m*Pr^(1/3)` factor cancels between the two. `jet_array_runner.py`
  calls `jet_array_impingement_nu` twice per point (at the point's own
  `Gc/Gj` and at `0`) and takes the ratio, rather than reimplementing the
  formula to compute an "expected" value.

  242 points score at bias +3.7%, RMS 9.6% against the implementation --
  close to the fit's own stated 5.6% standard error, given this is raw
  scatter rather than the fit's own training residual. One panel is a
  genuine outlier: `B(5,4)I`'s `z/d=1` series sits at the correlation's
  validity corner (`xn/d`, `yn/d`, `z/d` simultaneously at their lower
  bounds of 5, 4, 1) and underpredicts measured crossflow degradation
  there, error growing from +14% to +42% with `Gc/Gj` -- calibration ticks
  on that same panel read back within 0.9% of round numbers, so this reads
  as a real weak spot of the fit at its own domain corner, not a
  digitisation artifact.

  Axis calibration checked two ways per panel (tick positions against
  round numbers, corners against the nominal box) rather than relying on
  `verify.py`'s frame-slope check, which assumes a physically-drawn
  distortion line spanning a real data range and throws on the
  near-origin coordinates a tick-mark click produces.

  Figure 6 covers the inline pattern only, per its own caption --
  `florschuetz_1981_staggered()` remains unvalidated against scatter data.
  The paper's Figure 7 gives a staggered/inline Nu ratio for a handful of
  geometries, which could support an indirect check later, but is not
  digitised.

- **`han_park_1988_angled()`, a third rib correlation set** for angled
  ribs in broad-aspect-ratio rectangular ducts (Han and Park 1988, Eq.
  4.17/4.18): `alpha` 30-90 deg, `W/H` 1-4. Han's original power-law
  schema couldn't express this equation's two shapes -- `R` is a
  quadratic in `alpha`, and `G`'s `alpha`/`p_e` exponents switch on
  whether the channel is square -- so `RibCorrelationSet` gained two new,
  opt-in discriminators (`RAlphaShape`, `GShapeModel`). Both default to
  the existing plain power-law behaviour: `han_1988_orthogonal`'s and
  `rallabandi_2009_high_re`'s full pre-existing test suites pass
  unmodified, byte for byte.

  Both switches are genuine discontinuities the source states, not
  numerical artefacts -- up to 62.5% in `R` at `alpha=90`, `W/H=4`; up to
  27% in `G` at `W/H=1`, low `alpha`. Neither is smoothed, since no smooth
  transition is stated anywhere in the confirmed extraction; a solver
  that traverses either boundary exactly needs its own guard.

  Scored against its own source, figure 4.47: `R` at RMS 10.5% over 39
  points, `G` at RMS 8.8% over 115 (representative geometry, the cloud
  carries no legend). Independently cross-checked against figure 4.51's
  parallel-rib classes -- a different paper, Han et al. (1991) -- at RMS
  8.5%. The other seven rib shapes on that figure (crossed, V, lambda)
  are deliberately not scored: this equation has no way to represent rib
  shape beyond angle, and scoring a shape mismatch would repeat the
  mistake extraction item 23 already recorded once.

  Along the way, the validation harness gained a second scoring path (`R`
  indexed by `alpha` rather than `G` bisected to a target `e+`, since `R`
  carries no Reynolds-number dependence in every set implemented so far)
  and corrected two more stale `scores` fields on figure 4.51's parallel
  classes, silently refused since before this set existed to score them.

- **`rallabandi_2009_high_re()`, a second rib correlation set** for 45 deg
  sharp-edged ribs at Reynolds numbers an order of magnitude above
  `han_1988_orthogonal`'s range (30,000-400,000 vs 10,000-60,000). Named
  and user-selected, same as the existing set -- there is no
  auto-switching, matching #334's decision that users pick and define
  correlations explicitly, with provenance and validity visible.

  Scored against the 38-point digitised cloud from figure 4.193c: mean
  bias -2.8%, RMS 7.0%, independently reproducing the extraction
  document's own comparison (0.994 mean, 6.9% RMS) through the full
  `Re -> R -> f -> e+ -> G` chain rather than a direct `G(e+)` evaluation.

  Round-edged ribs at this Reynolds range are reported by the source to
  instead follow `han_1988_orthogonal`'s correlation, but the source
  calls that agreement "coincidental" rather than physically grounded, so
  it is not wired as an automatic fallback -- a caller with round-edged
  ribs chooses the correlation directly, the same as any other set
  choice.

- **Figure 4.48 digitised**: narrow-aspect-ratio channels (Eq. 4.19), both
  panels, plus a reference curve labelled `[Ref. 3]` and its error bars --
  16 series, 111 points.

  **Extraction item 41 resolved.** The `W/H=2/4` correlation, at Eq. 4.19's
  branch boundary, free-fits `(e+)^0.3302` -- matching the wide branch
  (`n=0.35`, 6-10% off) far better than the narrow one (`n=0.258`, 22-69%
  off). Han's own figure uses the wide-branch form at the boundary; the
  discontinuity itself remains, but which side it falls on no longer does.

  **A cross-paper gap surfaced, and nothing is scored across it.** The
  `[Ref. 3]` curve matches Eq. 4.17 and Eq. 4.18 (already confirmed on
  figures 4.46/4.47) to 4.16% and 3.70% RMS -- a third independent
  confirmation of both. But `Eq. 4.18`'s form (`C=2.24`, `n=0.35`)
  disagrees with `han_1988_orthogonal`'s own equation (`n=0.28`) by 11-21%
  at the same operating points, a real difference between two published
  papers rather than a model defect. Every series in this figure is
  therefore unscored, for three independent reasons recorded in the
  metadata: the new data falls outside the implemented set's valid `W/H`
  range, the reference curve is a drawn line, and scoring the overlap
  would report the cross-paper gap as model error.

  **A misread axis, corrected before filing.** Panel (a)'s `R` axis was
  initially logged as log; it is linear. The tell was the digitised
  frame's own picks -- `1.5, 2.5, 3.5, 4.5` are half-integer minor ticks,
  which a log axis cannot produce.

  **A third confirmation of Han's ~6% accuracy band**, this time from the
  source's own printed error bars: four reconstructed centers give
  half-widths of 6.5-8.2%, mean 7.2%. One bar's upper end was initially
  obscured; the rescan agreed with the earlier partial reading to 0.7%.

- **Everything is namespaced by source.** A figure number is not a unique
  identifier -- two books can each have a figure 4.46 -- so `cards/` and
  `plots/` are now keyed by source as `data/` already was:

  ```
  data/han2012/fig4.46_R_eD0.047_pe10_wh1.csv
  cards/han2012/fig4.46_R_scatter.yaml
  plots/han2012/fig4.46_redrawn.png
  ```

  `--figure` takes a bare number while it is unambiguous and the
  qualified `source/figure` form always. A bare number matching more than
  one source is an error naming the alternatives, never a silent pick.

- **Redrawn figures are committed.** All six live in
  `validation/cooling/plots/<source>/`, so a diff shows the plot rather than
  requiring a reviewer to regenerate it, and a change that moves points
  shows up as a changed image. Each carries its source citation in a
  footer naming the book and the primary paper behind the figure.

  `--all` redraws the whole dataset in one command, and a test pins that
  every figure has a committed plot -- the way this goes stale is filing a
  figure and forgetting to regenerate.

- **Figure 4.47 digitised**: the angled-rib correlation of Han and Park
  (1988), both panels plus four frames, 186 points.

  Equation 4.18 confirmed -- the pooled cloud free-fits `(e+)^0.3453`
  against the printed exponent of **0.35**. Equation 4.17 confirmed too:
  the quadratic in rib angle is reproduced with a minimum near 60-68
  degrees, mean bias -3.0% against Han's stated 6% band.

  **Pooling is forced here, not chosen.** Both panels draw several
  distinct symbols and carry no legend at all, so there is nothing to
  attribute them against -- unlike figures 4.46 and 4.51, where a legend
  made attribution hard but possible. It costs little: the left panel's
  ordinate is already normalised by `(P/e/10)^0.35 (W/H)^m`, and on the
  right only the exponent is comparable anyway.

  This figure is the source of extraction item 10 -- `G_bar = 1.2 G` is
  printed on its dashed line -- which figure 4.51 later measured at 1.2193.

  Not scored by any implemented set: `han_1988_orthogonal` binds 90 degrees
  only, and this IS the angled correlation. It is the acceptance data for
  a `han_park_1988_angled` set (#334).

- **Figure 4.51 digitised, both panels**: nine rib configurations in `G`
  and `G_bar`, two drawn reference lines, and four panel frames -- 24
  series, 135 points.

  **Extraction item 10 is now confirmed against data.** `G_bar/G` measured
  over the seven undisputed classes gives **1.2193** at 90 degrees against
  the printed `4.5/3.7 = 1.2162`, agreeing to 0.3%. Until now that ratio
  rested on two printed coefficients and a figure label. The angled
  configurations average 1.155, 4.5% below, against 2.6% scatter -- which
  is what justifies the runner refusing to convert `G_bar` off 90 degrees.

  **Two classes are unresolved and marked disputed.** `45 deg //` and
  `45 deg x` return `G_bar/G` below 1, which is impossible, and their two
  panels hold what appear to be the same marks. Three readings survive and
  the digitised data separates none of them; resolving it needs the
  primary paper rather than the reprint.

  Frames now declare `frame_orientation`. A vertical frame must be fitted
  across the other axis, and it measures shear -- which a horizontal frame
  cannot see. The two y-axis frames also pinned this figure's panels as
  aligned to 0.53%, which is what validated transferring the abscissa to
  the upper panel, whose x ticks are unlabelled.

- **Digitised figures can be redrawn from the committed data.** `uv run
  python -m validation.cooling.plot --figure 4.46` puts the data back on
  axes so it can be held next to the scan -- the check numbers cannot
  make, since a self-consistent series can still have the wrong shape or
  be swapped with another.

  It reads the same metadata the verifier and scorecard use, so a plot
  cannot drift from what the checks believe. Series now declare
  `x_axis_type`/`y_axis_type` and plot limits in their figure card, and a
  `figure`/`panel` assignment; a test pins that every series has them,
  since defaulting a linear figure to log looks plausible and is wrong.

  The equation a figure PRINTS is overlaid on the line it DRAWS, making
  that gap visible directly. Disputed class labels are marked in the
  legend.

- **A disputed class label is now a first-class state, not a comment.**
  `class_confidence` records how much a digitised series' SYMBOL-CLASS
  label is trusted -- `confirmed`, `provisional` or `disputed` -- as
  distinct from `confidence`, which is about the coordinates.

  A disputed series still pools, because its coordinates are sound and
  nothing in the current scoring depends on which symbol produced them.
  The runner refuses it the moment a correlation set actually binds
  geometry, since then the label IS an input. The verifier reports it on
  every run, and a test asserts the marking, the report, and the presence
  of a specific question to put to the source page. A fact that lives only
  in a metadata comment has a half-life; one that prints on every
  verification does not.

  One series is currently disputed: figure 4.46's `eD0.047/pe10/wh2` upper
  panel marks pair with their own R counterpart on only 2 of 4. The
  digitised data cannot settle it -- identification by elimination put the
  best alternative at 1 of 4 -- so it carries the question to ask of
  p. 376 rather than a guess.

- **Figure 4.46 digitised: `han_1988_orthogonal` validated against its own
  source.** 23 series -- both panels, all ten symbol classes, the three
  drawn correlation lines -- plus both panel frames as distortion
  standards. Pooled, the lower panel sits at **+0.1% bias across 64
  points** against the printed `R/(P/e/10)^0.35 = 3.2`.

  **Scoring is pooled, not per-class.** Symbols on the figure overlap and
  a mark cannot always be assigned to its class; some runs appear in only
  one panel. That costs nothing here, because the set carries zero
  geometry exponents in both `R` and `G` -- Han normalised `P/e` out of the
  R ordinate himself and `W/H` drops at 90 degrees -- so every class must
  land on one curve. Class labels are kept as provenance, marked
  provisional, rather than as something the scoring leans on.

  **Han's stated accuracy is a standard deviation, not a 95% bound.** The
  book says "within 6% for 95% of the data"; the digitised cloud puts 70%
  inside 6% with sd 5.9%, which is 1 sigma. Harness tolerances are set
  from the measured figure. A test pins it, because the alternative
  reading would justify a band twice as wide, and a band twice as wide is
  how a real error hides.

- **Digitised series are cross-verified against a figure card.** `uv run
  python -m validation.cooling.verify` checks every series against what the
  printed axes and equations say, read off the page independently of where
  the digitiser put the points. The two are separate channels and the check
  is their disagreement; which one is wrong is raised for a human rather
  than decided by the tool.

  It catches a dropped axis multiplier (the real figure 4.54 defect, two
  decades), a mis-calibrated span, two curves on one figure swapped, a
  drawn line that fails to reproduce its own printed equation, and
  double-picked or missed marks. Each was verified by injecting the bug.

  A bound that cannot be read off the scan is recorded as null and left
  unchecked rather than invented -- figure 4.193c's ordinate continues
  below its lowest labelled tick. The verifier's first run caught that
  mistake in the card itself.

- **A scored validation harness for the cooling correlations**, on the
  pattern `validation/junction/` uses. `uv run python -m
  validation.cooling.scorecard` reports MAE, RMSE, bias and in-band counts
  per digitised series.

  It scores THROUGH the model rather than around it: the runner bisects the
  Reynolds number until `evaluate_rib`'s own `e+` reaches the digitised
  abscissa, so friction factor, `e+` and `G` are all exercised. Recomputing
  the formula in the harness would score the model against a second copy of
  itself, which is the failure #333 exists to catch.

  The first source is eight series digitised from Han, Dutta & Ekkad
  (2012), figures 4.53, 4.54 and 4.193c. `han_1988_orthogonal` meets its
  own correlation line at 3.5% RMSE and the measured points at 5.4%,
  against Han's stated 6%. The digitised coordinates are committed; the
  scans stay gitignored.

  **Refusing to score is itself a result.** A set outside its numeric
  validity is still answering, and is reported as extrapolating. A set
  asked about a configuration it has no binding for is not:
  `han_1988_orthogonal` is 90 degree only, and scoring it against 45 degree
  rib data produced a 26.7% bias that reads as model error when it is the
  wrong correlation entirely. The runner declines and prints why.

- **Ribbed surfaces are selectable in the GUI again.** The node type returns
  with the fields the rebuilt correlation needs: rib geometry, the channel
  aspect ratio `W/H`, `n_ribbed_walls` as a 1/2/4 choice, and the
  smooth-wall Nu multiplier with an inline note on what it is for.

  The mapping tolerates an older saved network: a rib node stored before the
  newer fields existed loads with them defaulted rather than raising an
  `AttributeError` a user cannot act on.

  Dimpled, pin-fin and impingement stay unavailable and still reject on solve.
  Their parameter blocks remain in the inspector, unreachable from the
  dropdown, so a pre-0.7.0 network's parameters stay visible rather than being
  silently dropped on the next save.

- **Ribbed channels are back**, built on the parametrised correlation sets.
  `RibbedModel` carries a `RibCorrelationSet`, the rib geometry, the channel
  aspect ratio, and `n_ribbed_walls`.

  **The correlation gives the ribbed side; the element does the combining.**
  How many walls carry ribs is a design choice rather than a property of the
  correlation -- the channel is the internal wall by definition, and the other
  side of it is a different channel. `n_ribbed_walls` of 1, 2 or 4 sets the
  area weighting, with 2 meaning two OPPOSITE walls, the configuration Han
  measured. Three is rejected rather than interpolated: it has no unambiguous
  geometry.

  **Friction needs no wall weighting.** The `f` the roughness function is
  defined against is already the equivalent four-sided channel value, so the
  correlation returns what the channel needs and `ChannelElement` uses it
  directly. Nothing multiplies pipe friction -- correlations own their `f`,
  which is what #331 established after a round-trip through a restated
  friction factor moved a drop by -24.7%.

  **A documented gap, with the knob that closes it.** Using the plain smooth
  correlation for the smooth walls gives `h_s/h_r` of about **0.42**, where
  Han's own reported channel average implies **0.70** -- because ribs enhance
  the adjacent smooth wall by 10-50% as well, which no correlation here
  covers. The channel average is therefore about **20% below** Han's
  measurement for the two-ribbed-wall square case.

  That is recorded rather than papered over with an invented constant.
  `smooth_wall_Nu_multiplier` is the supported way to close it: about **1.67**
  reproduces Han. The size of the gap and the effect of the knob are both
  pinned as regression tests.

  The Jacobian is simple because `R` carries no `e+` term, so `f` does not
  depend on Reynolds number and `df/d(mdot)` is zero. The drop is odd in mass
  flow and its derivative therefore even -- it does not flip sign at zero,
  which is the direction that puts a sign error exactly at the crossing if
  reversed. Verified against central differences to 1e-10, including reversed
  flow.

- **Parametrised rib correlations.** Rib correlations are now data rather than
  code: a `RibCorrelationSet` carries the coefficients, each term's
  **normaliser**, the advisory validity band, the stated accuracy, and the
  provenance of all of it. `han_1988_orthogonal()` ships Han (1988) for 90 deg
  orthogonal ribs, taken from the confirmed extraction in
  `validation/cooling/extractions/han_ribbed.md`, and reproduces its worked
  values -- `R = 3.2000`, `f = 0.04576`, `e+ = 71.1`, `G = 12.21`,
  `St_r = 0.00968`.

  **Users can supply their own**, which is the expected path for real
  hardware: no published correlation is precise enough for a specific rig. A
  set carries a `RibProvenance` of `Extracted`, `Fitted` or `User`, so a
  coefficient tuned on a rig is structurally distinguishable from one
  transcribed out of a paper rather than distinguishable only by convention.

  Three design points, each earned by a defect in what preceded it:

  - **The normaliser is data, not a convention.** `3.2 (p/e/10)^0.35` and
    `1.4294 (p/e)^0.35` are the same function; applying one source's constant
    under another's convention is wrong by `10^0.35 = 2.24x` **uniformly**,
    which never looks like a trend and so survives review. A set whose
    reference is missing fails to load rather than defaulting to 1.0.
  - **Validity is advisory.** `RibResult.extrapolated` reports; nothing
    refuses. A band belongs to the source's rig, not the caller's hardware.
    Malformed *parameters*, by contrast, are rejected loudly by
    `validate_rib_set` -- they are mistakes, not operating points.
  - **`evaluate_rib` never throws.** Reverse flow, zero flow and extreme values
    are guarded smoothly rather than clamped, because the solver probes states
    that are not physical and a hard `abs()` or `max()` puts a kink in the
    Jacobian exactly where Newton iterates. The `e+` floor distorts the
    physical region by 1.4e-07 at `e+ = 1000`.

  `R` carries no `e+` term by construction, so **`f` is independent of
  Reynolds number** for this family and `df/d(mdot)` is identically zero. True
  of both correlations in scope. A correlation whose `R` varies with `e+`
  cannot be expressed here and would need its own implicit solve.

  Not yet wired into `ChannelElement`: that needs the ribbed-versus-smooth wall
  weighting, and is tracked in #334.

- **`ChannelResult.ddP_dvelocity`** [Pa*s/m], for chaining a channel
  correlation's pressure sensitivity onto a caller's own mass flow.
  `ddP_dmdot` is taken w.r.t. the mass flow through the correlation's internal
  flow area -- the pin-array minimum section, the jet holes -- so chaining it
  directly is wrong by the area ratio, **45x** for a 25 mm channel over a 3 mm
  pin array. The impingement routine is driven by per-jet mass flow rather
  than a velocity and so leaves the new field unset; `ConvectiveSurface`
  carries the split that converts it.
- **Pressure boundaries declare how they couple to a duct.**
  `PressureBoundary(coupling=...)` takes `"auto"` (default), `"static"` or
  `"total"`. A duct meeting a boundary either keeps its exit dynamic head (full
  recovery, a diffuser) or loses it to the plenum (a bare exit discharging to
  atmosphere) -- the previous behaviour was always the first, which is an ideal
  diffuser on every outlet.

  Pinning stagnation pressure at an outflow caps the mass flux at the sonic
  value FOR THAT PRESSURE, a constraint the physical problem never imposed.
  Measured on a 1.04 kg/s combustor duct at Tt = 1738 K: under total coupling
  D = 0.10 m demands 1.34x an impossible ceiling and will not solve; under
  static it runs at M = 0.752 and solves at Pt_in = 155 kPa. A genuinely
  oversized demand (D = 0.08 m) still fails either way, so real limits surface.

  `auto` infers from flow direction -- inflow takes total, outflow takes
  static. It cannot decide for a boundary whose flow reverses, or that serves
  both roles; that is what the explicit setting is for. Whoever builds such a
  network knows whether the exit is a plain opening or a diffuser, and the
  solver does not.



- **The momentum chamber closed on incompressible Bernoulli.**
  `MomentumChamberNode` enforced `Pt = P + 0.5*rho*v^2`, which under-predicts
  the stagnation rise by 11% at M = 1.0 and 34% at M = 1.4 -- and does so
  without bound past sonic, so a station driven transonic got a closure that
  was both wrong and unable to say so. It now closes on the isentropic
  stagnation state from entropy conservation with variable cp, which
  `stagnation.h` already provided.

  It is also a large convergence win where a root exists: the combustor network
  behind #351 solves in 1.2 s where it took 18.3 s, because the correct
  relation is 1.8 - 2.4x steeper through the transonic region and pushes Newton
  out of that territory on physics rather than on a tuned penalty.

  **The junction reference head follows.** `K` is tabulated as `dPt/q_dyn`, so
  the reference has to be whatever the ports call their dynamic head. While the
  node was incompressible the two agreed and their errors cancelled; correcting
  the node alone would have left the two halves of a junction disagreeing about
  one quantity. Both junction models now take the common port's own isentropic
  rise, computed from its `P`, `rho` and `u` rather than from `Pt - P` -- the
  same value at the solution, but it keeps the dependency set unchanged and so
  preserves the one-`Pt`-per-row structure the assembled Jacobian relies on.
  Taking it from `Pt - P` instead measured 7x slower.

  `MpceGeometry` gains `gamma` per port, following the kernel's own rule that
  the equation of state is evaluated at the call site where the mixture is
  known. For an ideal gas `a^2 = g*P/rho`, so `M^2 = rho*u^2/(g*P)` and no
  temperature is needed. Left at zero it falls back to the incompressible head,
  so a caller that cannot supply one keeps its previous behaviour.

  `MomentumChamberResult` gains `d_res_dT`: the incompressible form touched T
  only through `rho`, weakly enough that the assembled Jacobian tolerated its
  absence, which the stagnation closure does not. All derivatives are
  FD-verified to ~1e-9, including reverse flow and the `m_dot = 0` symmetry
  point.

  Nine tests asserted `0.5*rho*v^2` directly. Each now asserts the compressible
  rise **and** that it exceeds the incompressible value, so the relationship
  between the two is pinned rather than just the new number. Three others had
  lost their subject entirely -- the closure removed the pathology they were
  built around -- and were re-derived rather than relaxed:

  - the barrier-weight control now measures over a population, because the
    traced case stopped discriminating at any size across twelve decades;
    re-measured, the derived weight still beats the fixed fallback on 2 of 100
    scaled junctions, so the machinery stays.
  - the barrier fixed point is pinned at the element level, because the traced
    case now converges at every alpha from 1e4 to 1e9 and never goes
    wrong-direction, leaving no stalled iterate to read it from.
  - the assembled-Jacobian check now gates on whether each finite difference
    has converged, instead of an absolute floor that made it test roundoff on a
    structurally zero term.

  See issue #357.


- **The Fanno march handed Newton a staircase, and the choke barrier hid it
  behind a bigger one.** Two defects in the compressible channel's drop, both
  in the region a solver has to cross.

  When the flow chokes before the duct end the march truncates mid-step.
  `L_choke` was interpolated within the breaking step -- a comment above the
  loop explains why -- but the **outlet state was not**, so it stayed snapped
  to the march grid. `channel_compressible_mdot_and_jacobian` builds its drop
  from that outlet and differentiates it by central finite differences with a
  1e-6 relative step, about 0.1 Pa against risers of ~55 Pa. Scanning `m_dot`
  through the choked band, the gradient alternated sample to sample by 22%,
  and the reported Jacobian disagreed with a coarse finite difference by
  7 - 62% inside the band while matching it exactly outside. The outlet is now
  interpolated to the same station as `L_choke`.

  The barrier that keeps Newton out of infeasible flow then turned out to be
  compensating rather than correcting. The truncated drop covers the marched
  part only, so it collapses as choking moves upstream -- slope +13934 just
  below onset against -86257 just above, a sign reversal in the base term
  itself. The linear barrier did not remove that; it outshouted it, and the
  mismatch was the 74.5x derivative step at onset. The drop is now
  extrapolated to full length before the barrier is added, which is continuous
  through onset, and the barrier is quadratic in the truncated fraction so it
  switches on with zero slope. Monotonicity matters as much as smoothness
  here: a quadratic barrier alone leaves the drop decreasing with rising
  `m_dot`, which is exactly the spurious-root condition the barrier exists to
  prevent.

  Unchoked results are unchanged to the printed digits; RK4 was already
  grid-converged by `n_steps = 25` there. Pinned by
  `tests/test_fanno_choke_smoothness.cpp`, whose three cases were each
  falsified against the specific defect they cover. See issue #356.

- **Every compressible solver returned a state that could not be mixed.**
  `State::X` and `State::Y` are both public members, and only `set_X()` /
  `set_Y()` keep them in sync. All seventeen state-producing sites in
  `compressible.cpp` assigned the `X` member directly, so `fanno_channel`,
  `fanno_channel_rough`, `nozzle_flow`, `nozzle_quasi1d` and
  `solve_A_eff_from_mdot` handed back states with a populated `X` and an
  empty `Y`.

  Nothing caught it because every `State` property getter -- `h`, `cp`, `rho`,
  `mw`, `mu`, `gamma`, all of them -- reads `X`. Such a state computes every
  property correctly and looks healthy. `Y` is read in exactly one place a
  caller is likely to reach, `mix()`, which indexed it unchecked: passing a
  solver's own output straight into `mix()` segfaulted, with no misuse of the
  API anywhere.

  The producers now build through `set_X()`. `mix()` additionally checks the
  two sizes agree and throws naming `set_X()`, so a state built by hand the
  same way fails with something actionable instead of a crash. The four
  state-assembling sites at the pybind boundary were converted too.

  One visible consequence: `set_X()` normalises, so a solver given an
  unnormalised `X` now returns a normalised one. Previously the returned `X`
  echoed the raw input while every property was computed from the normalised
  composition, so the state disagreed with itself; it now does not.

  Python was never exposed: the `X` property setter on the binding has always
  synced `Y`. `MixtureState` maintains its own invariant correctly, and
  `solver_interface`'s `Stream` carries mass fractions only, so neither was
  affected. Tracked as #352; making the invariant structurally unbreakable is
  still open there.

- **Every example runs again, and CI now notices when one stops.** The 0.7.0
  correlation removal broke five examples, and the 0.6.0 `pipe_* -> channel_*`
  rename had quietly broken three more that nobody had run since. All eight
  are fixed.

  `scripts/run_examples.py` already ran every example headless and exited
  non-zero on failure; it just was not wired to anything. It now gates
  `status-check`, so an example that stops running blocks the merge.

  **The job runs each example rather than importing it.** An example calling a
  removed function from inside a function imports cleanly and only fails when
  the function actually runs -- measured on four of them, an import check
  caught two.

  Two examples are beyond a rename and are listed in the runner's
  `KNOWN_BROKEN` map with their reason and issue (#351): one demonstrates a
  `CombustorNode(pressure_loss=...)` callback that no longer exists, the other
  builds a merging `MomentumChamberNode` topology that `FlowNetwork.validate`
  now rejects. Both are still run, so one that starts passing again fails the
  job asking for its entry to be removed, rather than sitting on the list
  forever.

- **Three C++ examples segfaulted and nothing noticed.**
  `stagnation_example`, `compressible_example` and `combustion_state_example`
  crashed on a null read inside `combaero::mix()`. They assigned the public
  `State::X` member directly instead of calling `set_X()`, so `Y` stayed
  empty and `mix()` indexed it. They compiled and linked cleanly throughout.

  The examples were built by CMake but never registered as tests, so building
  them was the only thing ever checked. Each is now an `add_test` smoke case,
  which puts them in `ctest` and so in CI across the existing matrix.

- **`_RibbedChannelResult` carries the wall-coupling derivatives.**
  `dh_dmdot`, `dh_dT`, `dT_aw_dmdot` and `dT_aw_dT` were missing, so a ribbed
  channel coupled to a wall fed the solver an incomplete Jacobian. Nineteen
  unit tests missed it because none of them coupled a wall; re-enabling the
  skipped coupled-liner example test is what surfaced it.

- **`ChannelResult::ddP_dvelocity` is retained** and still set by
  `channel_smooth`, though its only consumer went with the pin-fin path. It is
  correct, tested, and is the interface the re-added correlations will use.

- **`ChannelElement.diagnostics` reported a friction factor from a different
  correlation than the one driving the residual.** It hardcoded Haaland
  (rough) or Petukhov (smooth) regardless of the element's `friction_model`,
  so a channel set to `petukhov` reported **0.016834** while its own pressure
  drop used **0.012928** -- a 30% misreport. It now asks the same C++
  dispatcher the residual uses. The laminar branch stays `64/Re`, which is
  Poiseuille rather than a choice of correlation.

- **`dimple_friction_multiplier` no longer warns about a Reynolds range it does
  not have.** It shared a validator written for `dimple_nusselt_enhancement`,
  which really is a power law in Re, so it emitted "Re_Dh is outside validated
  range [10000, 80000]. Extrapolating power-law correlation" -- naming the
  Nusselt function, about an extrapolation that cannot occur. Geometry limits
  on `d_Dh` and `h_d` are still enforced, and the Nusselt path is unchanged.

- **A test that could not fail has been replaced.** `test_dimple_jacobians`
  compared the analytic derivative against a central difference of the same
  function; both were 0.0, so it held regardless of the implementation. The
  Nusselt half was real and is kept. The friction half is now a test that the
  multiplier does not vary with Re, which goes red if a Re dependence is added.

- **`ChannelElement`'s Jacobian was wrong for pin-fin and impingement
  surfaces.** Those surfaces vary their friction multiplier with mass flow,
  but the multiplier was handed to the C++ friction routine frozen, so the
  element reported `d(dP)/d(mdot)` short by **10.5%** and `d(dP)/dT` short by
  **13.1%** against central differences. Both now agree to 1e-10. Ribbed and
  dimpled were exact throughout -- their multipliers carry no Reynolds
  dependence, which is what localised the fault. Newton convergence on
  networks using either surface should improve; the converged solution for the
  default `friction_model="haaland"` is unchanged.

### Changed

- **An orifice plate between pipes charges its permanent loss, not its tapping differential (#645).** With a pipe on the side the flow comes from (which sets beta) and a wider pipe on the side it goes to, a metering (`ReaderHarrisGallagher`, `Stolz`) or `fixed` `OrificeElement` now loses ISO 5167-2:2003's permanent pressure loss, k = [s - C b^2] / [s + C b^2] of the tapping differential, s = sqrt(1 - b^4 (1 - C^2)).
  - **Before:** the whole tapping differential was lost, so the pressure recovery downstream of the vena contracta was missing: k is 0.90 at beta 0.3, 0.73 at beta 0.5 and 0.52 at beta 0.7, so the orifice charged up to about 2x its loss and pressure-driven networks under-predicted its flow.
  - **Unchanged:** an orifice discharging into a node with no pipe beyond it (a plenum, a boundary) still loses the tapping differential, and a plenum-fed orifice (beta = 0) is untouched. The discharge-hole correlations (Idelchik, McGreehan-Schotsch, Lichtarowicz) keep the old form in a pipe until their in-pipe forms land.
  - **The network's metering Cd is now evaluated at D-D/2 taps**, the taps ISO's loss relation is derived for, instead of flange taps: the flow moves by at most 0.2% up to beta 0.6 and 1.3% at beta 0.75.
  - **Compressible:** the nozzle expands to the tap pressure the recovered drop implies. ISO's relation is incompressible; `recovery_in_range` is False where p_tap/p1 < 0.75.

- **A `ThermalWall` passes the two streams' exchange, epsilon-NTU (#627).**
  `Q = eps(NTU, C_r) C_min (T_aw,a - T_aw,b)` at each stream's inlet,
  replacing the linear `U A dT`, which is its small-NTU limit.
  - **What the linear form did:** it moved more heat than the streams hold at
    NTU ~ 1, cooling a 1600 K stream to 298 K against a 400 K cold inlet (#551,
    part 1). It also kept a finite Q as one stream stopped, which sent 10 of
    40 stress hot-cold networks to a stalled hot side below 0 K. All of them
    now solve, with every outlet between the inlets.
  - **The heat shift:** the suite's 157 converged walls (NTU 0.01-1.06) move
    less heat: median 2.8%, 18% at the 90th percentile, at most 57% (the
    #551 case, where the linear form moved 2.3x the exchange).
  - **Capacity rates:** `C = |m| cp` with cp the secant over the inlet span
    (new `cp_secant_mass_and_jacobian`), so `C_min dT` is the exact maximum
    heat. A node side is a reservoir.
  - **New `ThermalWall.orientation`** (`"aligned"` default, or `"opposed"`):
    the geometry, from which the flow signs decide parallel or counter flow.
    The arrangement therefore follows a reversing flow. It is also passed
    through the GUI's thermal wall schema.
  - **New C++ kernel:** `wall_exchange_and_jacobian`, with analytic
    derivatives including `dQ_dC_a` and `dQ_dC_b`. `WallCouplingResult` gains
    `dQ_dC_a`, `dQ_dC_b`, `UA` and `effectiveness`.
  - **New solution keys:** `{wall}.UA`, `.NTU`, `.effectiveness` and
    `.counterflow`.
- **The LM fallback also follows a hybr "success" at a non-root** on any
  network, not only compressible, tee, multi-port or merge ones. An 8-row
  impingement array stopped at |F| = 1.7e-3, where LM from hybr's best
  iterate reaches 1.5e-10.

- **A converged root on which an element destroys entropy is now refused (#598).** An element whose second-law verdict is `"fail"` (`sigma_norm_high < -1e-8`: it creates flow work even on its most favourable heat path) demotes the solve to INCONSISTENT, `__success__` False, with the element named in the message. The check was report-only since #545. On the stress sweep the roots this newly refuses are an `EjectorElement` entraining more power than its primary jet carries, and an `IdealMomentumTeeElement` split with `loss_factor = 0` at Mach. On the incompressible closure the band is widened by that closure's own resolution, order M^4 (`sigma_norm_resolution`), so a split lossless in its own model is not refused for its ports' density spread. A ThermalWall's own check (`__second_law_walls__`) stays a report.

- Each node type now has one icon and one colour, shared by the palette and the canvas. 7 of 17 types used to differ between the two, and three pairs shared a glyph (mass/pressure boundary, effusion plate/liner, ejector/momentum chamber). Changed: the orifice is `tally-2`, the pressure boundary a gauge, the effusion liner four columns, the momentum chamber a circled dot (it no longer borrows the ejector's wind glyph), and the probe card shows its crosshair (#321).

- `ImpingementArray.summarize`'s `heat_transfer_outside_source` now flags a bypass more than 1 K off the supply only when it meets a staggered array or a row outside CR-3936's inline data. Inside that data, a hotter bypass is modelled rather than flagged. The GUI label reads "eta_r outside data" (#467).

- **Where the network ends is defined once, on the node (#545).**
  - **Before:** graph validation matched class names, the solver matched `PressureBoundary`/`MassFlowBoundary`, and the GUI export kept its own kind list.
  - **Now:** every node class says whether the network may end there (`is_boundary`: pressure, mass-flow and wall nodes) and whether mass and enthalpy cross the network's control volume there (`is_flow_boundary`: pressure and mass-flow nodes). Validation, the solver and the GUI export all read these attributes.
  - **Export:** the `is_boundary` column is now also true for wall nodes, which it previously left out.

- The species thermo now defaults to the C2 re-fit of NASA-9 (`thermo_model() == "nasa9_c2"`, #549). Outside the |T - T_join| < 0.05 T_join windows at the 1000 K and 6000 K interval joins, cp, h and s are identical to NASA-9. Inside them, cp moves by at most 4.8e-4 relative, h by 1.1e-5 cp T and s by 3.1e-5 R. dcp/dT no longer steps at the joins, so neither do the solver Jacobians that use it. `set_thermo_model("nasa9")` restores the published polynomials.

- **A wall-heated compressible channel carries its heat in the momentum (#551).** It used to be adiabatic Fanno flow, whose wall heat only raised the receiving node's Tt. That omitted the Rayleigh change: heating costs total pressure, cooling gains it.
  - **Now:** where the receiving node is a plenum or chamber fed by this duct alone (no wall of its own), the duct uses `fanno_rayleigh_channel_flow`, with that node's Tt as its exit Tt.
  - **Second law:** the four stress hot/cold pairs that failed it now pass.
  - **Into a node that mixes other streams, or a boundary:** the duct stays adiabatic and keeps the "Heating not in the momentum" flag (#660). The Inspector text now says why.
  - **Results change for walled compressible ducts:** heated ducts pass less flow, cooled ducts more. A strongly cooled duct can now choke at its inlet.
- **The duct heat-transfer state is smooth through choke.** `duct_htc_static_state` is the new state the wall's htc reads. It equals `duct_static_state` below 99.8% of the sonic flux, then bends C1 onto 99.9% of it (M below ~0.97).
  - **Why:** a duct choked at its inlet sat exactly where the isentropic state's slope in G is infinite, and its wall coupling could not converge.
  - **New fields:** `DuctStaticState` gains `G_eff` and its partials.
- **A smooth-wall channel whose heat-transfer area exceeds its friction-wetted area (4 A L / Dh) by over 10% is flagged** (`heat_area_ratio`, `heat_area_exceeds_wetted`, and an Inspector note).
  - **Why:** heat over surface the friction never sees breaks Reynolds analogy. In a compressible duct a cooling gain can then outgrow the friction, and at a fixed exit Tt no flow matches a pressure drop.
  - **Where it showed up:** the stress harness's impingement crossflow segments carried 3.9x (the plate's width on a round duct). They now use their wetted area, which changes that family's draws.
- **A cold solve of a network with heat-carrying compressible channels runs in two stages:** those channels adiabatic first, then the exact problem from that root. The plain cold chain runs only if that fails. An impingement crossflow chain went from 23 s to 4.7 s.

- **A newly dropped GUI channel inherits its diameter (#570).**
  - **Before:** the canvas seeded `D = 0.1` as a user value, so a fresh channel never inherited.
  - **Now:** D is unset. The channel inherits from a sized neighbour, and with nothing to inherit from it falls back to a marked 0.1 m default (`area_source = "default"`, amber in the Inspector). A freshly built chain of channels still solves.
  - New `ChannelElement(default_diameter=...)` opt-in, with `DEFAULT_CHANNEL_DIAMETER`. The Python API does not opt in, so a script with an unresolvable channel diameter is still refused (#262).

- **Breaking:** the flow regime now defaults to compressible (#585). `DEFAULT_REGIME`, which every element and `MomentumChamberNode` uses when no `regime` is given, and the GUI's global regime for a new network are now `"compressible"`; incompressible flow is its low-Mach limit. A script that relied on the old default solves the compressible forms now, so its results change wherever Mach number matters. Pass `regime="incompressible"` to keep the old behaviour. A saved GUI project keeps the regime it was saved with.

- **Viscosity and thermal conductivity are C2 in temperature** (#611). The Omega*(2,2) collision integral is a quintic Hermite in (ln T*, ln Omega*) with tabulated node curvatures (exact quadrature for Lennard-Jones). The cubic Hermite stepped d2mu/dT2 by 0.1-10% at every node. The Eucken cv_vib clamp is a C2 smooth positive part (`CV_VIB_SMOOTH_HALF_WIDTH` = 0.002 R); the hard clamp stepped dk/dT by up to 4% (H2 355 K, CH4 178 K, O2 189 K). Values move by at most 2.3e-5 in mu and 1e-4 in k from 200 to 4000 K (air: 1.7e-6 and 2.6e-5).

- **Cold solves start from a resistance-network solve.** `NetworkSolver` seeds node pressures and element flows from the network solved as quadratic resistances (dPt = K m|m|, mass conserved), whenever every element offers an estimate: orifices, channels, pressure losses with a shipped correlation, area changes and lossless connections. Other networks keep the old start. On the stress `random` family, 67 -> 113 of 160 draws converge, in about half the time; `solver.resistance_seed = False` restores the old start.

- A compressible `"default"` solve of a network containing an impingement array starts from the outlet-referenced incompressible warm start, falling back to the cold path, auto-retry and all, if that fails. Elements declare it with `prefers_warm_start`, and only `ImpingementCrossflowElement` does. The stress suite's impingement family went from 37/40 converged in 167 s to 40/40 in 33 s, every previously converged root unchanged; a 10-row array solves in 0.5 s instead of 1.4 s, and with a bypass (#467) in 0.5 s instead of 6-11 s (#706).

- **Area change closures follow Idelchik (1966) and carry Mach through their control volume (#731).** **Network results change.**
  - **Sudden expansion:** diagram 4-1. That is 26/Re below Re 8, a C1 interpolation of the printed table to Re 3500, and Borda-Carnot above. The old alpha = 1 + 2.73/(Re + 10)^0.4 had no source: -50% bias against the table below Re 3500, and +0.069 dynamic heads above it at Re 1e4, even with no area change.
  - **Conical diffuser:** diagram 5-2, exact at every printed cell. The old fit read 3x at full angles of 10 deg or less.
  - **Compressible:** an expansion is the compressible Borda-Carnot step times the table's ratio to it; a contraction is a vena contracta re-expanding by the same step, which chokes. `regime="incompressible"` gives the M -> 0 limit.
- **`AreaChangeElement` residual is a mass flow, `m - m_calc` (#731),** like the compressible channel and orifice.
  - **Flow:** it rises monotonically with the drive and saturates at the sonic flow, C1 through zero drive and through choke.
  - **Equal areas are refused.**
- **The GUI's hidden tee-port steps run the network's regime (#731).**

- **Sharp and conical contractions follow Idelchik (1966) (#731).** **Network results change below Re 1e4 and for conical nozzles.**
  - **Sharp contraction:** diagrams 3-10 (Karev) and 3-9. That is 27/Re below Re 8, a C1 interpolation of the printed table to Re 1e4 (every cell exact), and 0.5 (1 - n) above. It replaces a constant 0.5 (1 - n) at every Re.
  - **Conical contraction:** diagram 3-7, zeta' (1 - n), with zeta'(full angle, l/D) from diagram 3-6. It replaces the unsourced C_con = 0.04 + 0.46 sin(theta)^1.5.
  - **API:** `sharp_contraction_zeta(n)` becomes `sudden_contraction_zeta(Re, n)` (returns `(zeta, dzeta_dRe)`), and `conical_nozzle_zeta` takes `(full_angle_deg, length_over_D, n)`. Neither has been released yet.

- **One naming convention for result keys (#773).** `P` is a static pressure and `Pt` a total pressure on every element, a key means one thing on every element, and the ejector no longer puts the unit in its pressure keys. The renamed keys are the columns of the CSV export and of `NetworkRunner.to_dataframe()`, and the keys of the solver's result dict, so a notebook that reads one by name needs the new name. The physics is unchanged: every other value of the GUI's palette networks is bit-identical before and after.
  - **Junctions report `P_jct` (static) and `Pt_jct` (total), every model.** The Mynard and constant-K tees reported their junction TOTAL pressure as `P_jct`; it is now `Pt_jct`, and `P_jct` is the junction static, `Pt_jct` minus the common port's dynamic head (the common port's static at a converged solution). The ideal momentum tee keeps `P_jct` (its static `p_jn`), drops the duplicate `p_jn` and adds `Pt_jct`, `P_jct` plus the common port's dynamic head. The constant-K tee did not report `Pt_jct` before. The solver unknowns follow the same convention, so a variable and a result key of one name are one quantity: the Mynard and constant-K tees' unknown is `{id}.Pt_jct` (it was `{id}.P_jct`, holding the total); the ideal momentum tee and the ejector keep `{id}.P_jct`, a static. An `initial_guess` under the old `{id}.P_jct` key still seeds a Mynard or constant-K tee, with a `DeprecationWarning`, until 0.8.0; a saved GUI file's tee guess moves to `Pt_jct` silently on load.
  - **Ejector:** `outlet_pt_pa` -> `Pt_out`, `p_c_star_pa` -> `P_c_star`, `p_b0_pa` -> `P_b0`, `omega` -> `entrainment_ratio`, `omega_critical` -> `entrainment_ratio_critical`. `p_py_pa` is gone: it repeated `P_jct`, the mixing-plane static. The ejector reports no `Pt_jct`: its mixing plane carries two streams at one static pressure and different velocities, so it has no single total pressure. The internal ejector API (`_ejector_huang1999`, the C++ bindings, `EjectorEntrainmentResult.omega`) keeps Huang's `omega`.
  - **Orifice:** `dp_tap` -> `dP_tap`, `dp_permanent` -> `dP_permanent`, like `dP`, `dPt` and `dP_drive`.
  - **`area_ratio` is the area change's F1/F0 only.** The effusion plate's approach area over hole area is `approach_to_hole_area_ratio`; a junction port's connected-element area over port area is `port_<i>_connected_area_ratio` (`port_<i>_area_ratio` before), and the GUI tee's hidden port step reports `step_<port>_connected_area_ratio`.
  - **The shaft speed is `shaft_speed_rpm`,** not `omega_rpm` (omega is 2 pi f, not a speed in rpm): the `VortexElement` keyword, attribute and result key, the GUI vortex node and the GUI solver setting. `VortexElement(omega_rpm=...)` and `.omega_rpm` still work with a `DeprecationWarning` and will be removed in 0.8.0; a GUI network saved with `omega_rpm` still loads, in the browser and through the backend.

- **Changelog entries are fragments in `changelog.d/`, not edits to `CHANGELOG.md`.** Every PR used to write its entry at the top of `[Unreleased]`, so any two open PRs conflicted there.
  - **Now:** a change adds `changelog.d/<id>.<category>.md` (see `changelog.d/README.md`). `scripts/changelog.py release X.Y.Z` folds them into a dated block at release time, merging them with what is already under `[Unreleased]` one section per category. It also moves the compare links on and deletes the fragments.
  - **Validation:** `scripts/changelog.py check` runs as a pre-commit hook.

- The built-in combustor pressure losses (`ConstantFractionLoss`, `LinearThetaFractionLoss`, `ConstantHeadLoss`, `LinearThetaHeadLoss`) now evaluate in C++, through the new `cb.combustor_fraction_loss` and `cb.combustor_head_loss`, which return `CombustorLossResult` with every partial. `PressureLossElement` builds its Jacobian from them instead of differencing the correlation. Values are unchanged (they agree with the Python forms to round-off). A user-supplied loss callable keeps the differenced path (#718).

- GUI: a constant-K tee's K values now default to Bassett 2001's analytical pair (dividing K5/K6, joining K11/K12 angle-corrected) at a new **Design Branch Share** (default 0.5), the resolved arm areas and the branch angle, instead of the fixed 0.4 / 1.0. Those fixed values suit only an equal-area tee: with a branch a quarter of the main area the branch K is about 3.5 at an equal split and 7 at a 0.7 share. A typed K still overrides its own leg; clearing it returns to the computed value, which the Inspector shows greyed. Saved graphs with typed K values are unchanged; graphs that never set K now solve with the computed pair.
- `ConstantKTeeElement.bassett_design_K(flow_direction, q_branch, psi, theta_deg)` returns that pair.

- **Dependency floors raised to what is tested: `numpy>=1.26` and `scipy>=1.15`.** numpy 1.24/1.25 and scipy 1.10 have no wheels for Python 3.12, the only supported Python, so the old floors could not be installed. With scipy below 1.15, one junction stress case also stalls at a different point. A nightly `dependency-floors` job now installs every runtime dependency at its floor and runs the fast suite and the examples. PRs that touch a `pyproject.toml` run it too. Examples and tests use `scipy.integrate.trapezoid`, so they also run on numpy 1.26.

- A compressible `EjectorElement` now uses the `"huang1999"` operating anchors by default, in place of Kracik & Dvorak's lossless pair: a lower critical back pressure P_c* and a dead head P_b0 3-5% above the measured breakdowns, where the old one ran 62-274% above.
  - Results change for every compressible ejector: the entrainment plateau ends at a lower outlet pressure, and the droop falls to zero much sooner.
  - An outlet above the new dead head has no forward-flow root, so some networks that solved before now fail. The solver's message names the dead head.
  - Pass `anchor_closure="kracik_dvorak"` (GUI: "Back-Pressure Anchors") for the old pair.
  - Incompressible ejectors are unchanged. Their compressible warm-start proxy keeps the old pair's wider window, which converges more often and only seeds the solve.

- **Ejector mixing area A_3 and the outlet pipe inherit from each other (#559).** The mixing duct and the pipe downstream are the same duct, so whichever side is set wins: a pipe with no diameter takes A_3 from the ejector, and an ejector with `mixing_area=None` takes it from the pipe. The 8e-4 m² default applies only when neither is set. The inlet faces never inherit (the nozzle exit and the suction annulus are internal geometry, not the feed pipes' areas). An inherited A_3 that is not larger than the nozzle exit is refused with a reason. New ejectors in the GUI start with A_3 unset; saved networks keep their stored value.
- The ejector's Inspector shows throat, nozzle exit and mixing duct as D | Dh | A blocks, like every other element.
- Fixes found in live testing: a cleared ejector throat or nozzle-exit area now means its default everywhere (it was shown as the default but validated as 0, blocking Solve); the computed counterpart of an entered D or A is shown as a rounded placeholder instead of a full-precision value that read as typed input; a field that cannot inherit no longer tells you to "untick Inherit".

- **Every network element declares the units of its residual rows; there is no default.** Each element sets `residual_scale_kind` ("p" or "mdot"), and junctions set `row_scale_kinds()`. `FlowNetwork.validate()` raises `TypeError` unless the declaration sits on the class that defines `residuals` or a subclass of it. A custom element, or a subclass that rewrites `residuals`, must declare. The silent "mdot" default had misfiled pressure rows three times (#614). Solver behaviour is unchanged: the suite's 26 124 solves take the same evaluations.

- **The nozzle inverse solvers and `fanno_max_length` refuse invalid input.** `solve_A_eff_from_mdot`, `solve_P_back_from_mdot` and `solve_P0_from_mdot` clamped a non-positive T0 to 10 K and a pressure to 1 Pa and returned an answer for that state (`solve_P0_from_mdot` with T0 = -600 K returned 1e18 Pa; `solve_A_eff_from_mdot` returned a negative area for a negative flow). They now raise `ValueError` (`std::invalid_argument`) for non-positive temperatures, pressures, areas or target flows. `fanno_max_length` raises for non-positive `T_in`, `P_in`, `u_in`, `D` or `f`, where it used to return a value without complaint. `nozzle_flow` itself still clamps, since network Newton iterates call it.

- **Every element shows its cross-section the same way: D | Dh | A (#559).** Channel, area change (per face), tee (per arm), discrete loss, combustor and momentum chamber now share one geometry block.
  - D and A are two views of one value: enter whichever you have and the other is computed (A = πD²/4; for a non-circular duct D is the area-equivalent diameter). Dh = D unless overridden on elements that store a Dh.
  - Inherited geometry is locked by default and shown greyed with its source. Untick **Inherit** to opt out: the current inherited value becomes your own and the fields unlock. Ticking it again clears the override. This replaces "edit to override" / "Reset to inherited", under which a seeded value silently counted as user-set.
  - An element with nothing to inherit is editable directly, with an amber hint; a nominal default (chamber 0.1 m², tee `F_C` 0.01 m²) is marked as such.

- **One regime default for scripts and the GUI (#585).** Every element constructor that takes a `regime`, the GUI backend's `SolverSettings.global_regime` and the GUI store now share `combaero.network.components.DEFAULT_REGIME`, and a test keeps them equal.
  - It is **incompressible** for now: with compressible as the default, 60 tests failed to converge.
  - **GUI:** a newly drawn network starts incompressible; before, the store defaulted to compressible. Saved networks keep their saved setting, and Solver Settings switches it.
  - The backend fell back to incompressible already, so headless runs are unchanged.
  - The switch to compressible is tracked in #585, with its acceptance test.

- GUI CSV export / `to_dataframe()` unit labels are ASCII: `rho [kg/m^3]`, `mu [Pa*s]`, `htc [W/m^2/K]`. They were `kg/m³`, `Pa·s` and `W/m²/K`, which Excel showed garbled because the file has no byte-order mark. Code that selects these columns by their full header must use the new spelling. The GUI still displays the typeset units. The unused `Q_fraction` (labelled W for a fraction) and `ratio` unit entries are gone.

- **GUI polish (#321):** the lossless connection's port labels are `i` / `o` (matching the channel and the palette legend) instead of the clipped `in` / `out`, and its title follows the channel card's typography. `scripts/run-gui.sh` binds the dev server to `127.0.0.1` instead of every interface.

- GUI: a thermal wall whose box was moved now draws its arrow through the box (node -> box -> node) as one smooth curve, so the box stays visibly attached to its wall. An unmoved wall keeps React Flow's own curve.

- **Species thermo below 200 K is an ideal gas down to 0 K.** Below the
  NASA-9 data, cp/R now fades over the whole [0, 200 K] (C2 at 200 K) to the
  rigid-rotor value, 5/2 for an atom (and H2, whose rotation freezes near
  85 K), 7/2 linear, 4 nonlinear, instead of levelling off at its 200 K value
  50 K below. h(0 K) is finite; `calc_T_from_h` resolves T down to a few K and
  throws below h(0 K). New constants `thermo::CP_R_RIGID_ATOM/_LINEAR/_NONLINEAR`;
  `thermo::TAIL_FRACTION_LOW` is removed.
- **No [200, 5000] K clamp on network node temperatures.** A node carries its
  T(h) as is: heated past 5000 K it is simply hot, energy conserved (the clamp
  discarded ~190 kW of 200 in a heated plenum). `__T_clamped__` now lists the
  nodes with flow at or below 0 K (more heat out than the stream holds as a
  gas), and such a node reports `{node}.no_gas_state` instead of properties.
  The mixer inverts its enthalpy exactly to 0 K (the 134 K linear
  continuation is gone) and continues through it on h's line.
- **A stagnant node takes up no heat.** Below `MIXER_HEAT_MDOT_FLOOR` (1 mg/s)
  the heat a node takes up fades C1 to none at zero flow
  (`mixer_heat_per_unit_flow`, new, with its derivative), so a dead end holds
  the temperature of what reaches it; the previous floor put a 100 W heater in
  a dead end at 1e5 K. A root inside the fade -- a heated node neither flowing
  nor stagnant -- is refused as INCONSISTENT and listed in `__heat_faded__`.

- **A second-law refusal of an ideal tee names its loss factor.** When `IdealMomentumTeeElement` runs with a `loss_factor` below about M²/4 at its port Mach, its low-Mach loss form creates flow work and the second-law check refuses the root (#598). The message used to give only a sigma value. It now says that the loss factor is the cause and that 1.0 is the parameter-free closure. The solver gains a generic `second_law_hint(sol)` element hook for this.

- **A failed solve says when a declared-direction junction is the reason.**
  `MultiPortChamberElement` and `ConstantKTeeElement` take a `flow_direction`
  and do not reverse a port. When a solve fails with such a junction holding
  a port at zero or against its declared direction, the message now names
  the junction and the port, and points at the boundary pressures,
  `flow_direction` or `IdealMomentumTeeElement`. Every one of the 8 stress
  `tee3` draws that stalled has a root with the ideal tee, with that port
  reversed; before, the message said only "not making good progress".

- `normalize_fractions` no longer floors every fraction at 1e-15. It clamps negatives to 0 and keeps absent species exactly 0. The network's state composition likewise clamps to [0, 1] instead of [1e-12, 1].
  - **Why:** the floors put a trace of all 15 species into every stream, and every mixture property then evaluated all of them. Dry air is 4 species.
  - **Speed-up:** the Fanno channel kernel runs 3.3x faster on air (13.5 to 4.1 ms per call). The slowest suite test went from 189 s to 65 s locally, and the full suite from about 245 s to 130 s.
  - **No numerical change:** no validation score moves at its printed precision. Mixture `cp`, `h` and `dcp_dT` skip absent species (bit-identical), and `isentropic_expansion_coefficient` evaluates the mixture cp once.

- **Orifice: the bore uses the shared D | Dh | A block, and β comes from the pipe the solver actually uses (#559).**
  - The 0.08 m default bore is shown greyed and marked as a default instead of being written into the field as if typed.
  - The pipe diameter behind β = d/D now comes from the solver's resolve pass, so an upstream pipe that *inherits* its D is picked up. Before, the Inspector saw only a D typed on the adjacent channel and showed no β otherwise. Its source is shown.
  - Warns when the bore is equal to the pipe D (no restriction; behaves as a continuous pipe) or larger than it (likely not intended); the solver skips the velocity-of-approach correction (E = 1) in both cases.

- **`scripts/run-gui.sh` rebuilds a stale C++ extension before launching.**
  The install is editable, so Python changes are live, but the compiled
  `combaero._core` is not. After a pull that touched C++, the backend died at
  import against the old build. The new `scripts/ensure-gui-build.sh`
  rebuilds when a C++ input is newer than the build, then checks that the
  backend imports cleanly, counting combaero's swallowed "Import failed"
  fallback as a failure. It stops with the error instead of launching a broken
  backend. `--check` reports without rebuilding.

- **Second-law entropy check tightened to 1e-8 (#545).**
  - The per-element internal-entropy tolerance `SIGMA_TOL` drops from 1e-6 to 1e-8, so the report now flags a manufactured total pressure down to the sub-ppm level -- 100x more sensitive than before.
  - Kept global, not per element: a sweep of 238 element instances (the suite's networks plus 50 random junctions) found nothing dipping below zero -- the smallest real loss was 1.5e-5 (an iterative junction) -- so no per-element floor was warranted, and 1e-8 keeps a ~10x margin over the ~1e-9 solver round-off.

- **The GUI backend and `NetworkRunner` default to a 30 s max solve time, as the GUI does.** A request without a timeout (an older saved file, a script, a direct API call) used to get 180 s; `NetworkRunner.solve` and `.sweep` defaulted to 180 s as well. The value is `gui.backend.schemas.DEFAULT_SOLVE_TIMEOUT`, and a test keeps the frontend store equal to it. `None` still means no limit.

- **A numerical exception at a solve's starting point no longer aborts `solve()`** (#552 follow-up). The attempt ends with outcome `ERROR`, a `RuntimeWarning` names the node or element, and the auto-retries try their other starts. Code bugs (`TypeError`, `AttributeError`, `KeyError`, `IndexError`) still raise. `solver.raise_internal_errors = True`, or `COMBAERO_SOLVER_STRICT=1`, restores raising; the test suite runs strict.
- **`__caught_exceptions__` / `__caught_exception_count__`** in the solve result record every exception the solve caught: the starting point and the trial points the root finder rejected, with phase, node or element (or the non-finite Jacobian rows and columns), type, message and innermost frames. This is a discovery log for stress runs.
- **`WallLayer` raises on an unknown material name.** It used to be swallowed at every evaluation, silently keeping the manual conductivity.

- `CombustorNode`'s unburned temperature (the `T_unburned` of theta = T_burned/T_unburned - 1 that a theta-sourced `PressureLossElement` reads, and that loss's reference temperature) is now the enthalpy-mixed reactant temperature, the same mix the C++ combustor context defines theta with. It was a mass-weighted Tt. With air at 700 K and CH4 at 300 K (phi about 0.86) it reads 653.3 K instead of 681.0 K, which moves theta by about 4% (#718). Its Jacobian now carries each inflow's composition as well as its flow and Tt, all from the C++ mixer's analytic partials.

- **`NetworkSolver` assembles the Jacobian only when the root finder asks for it (#489).**
  - **Before:** the residual callback returned (F, J), so every evaluation built the state-sensitivity relay (dT/dx, dY/dx) and the sparse Jacobian. hybr uses J only at the start and on restarts: 74 of 401 evaluations on a 36-case mixing study.
  - **Now:** F-only evaluations skip the relay. A Jacobian request at the point F was just evaluated reuses that evaluation's element derivatives and rebuilds only the relay, so the answer is the eager Jacobian and asking costs no second evaluation.
  - **Results unchanged:** on 292 networks (mixing plenums with compressible ducts, hot impingement arrays, 240 random junctions) the solutions are bitwise identical and the evaluation counts are the same.
  - **Faster:** impingement arrays 2.1x (3.3 -> 1.6 s), random junctions 24% (7.2 -> 5.5 s). Mixing networks are unchanged (2.3 s); their cost is the Fanno kernel.
  - **Penalty reference:** a rejected probe's penalty still references the residual at the root finder's accepted point, which is the point it asks J at. Referencing the last trial point instead turned a converging MPCE cold start into a warm-start retry.

- **`ImpingementPlateElement` refuses wiring it does not model (#465).** The
  Florschuetz configuration is spent air leaving down the crossflow channel.
  The plate now refuses three other configurations at set-up:
  - a bypass or initial crossflow feeding the chain (#467);
  - a gap draining through the target, as in impingement-effusion or film
    (#468);
  - two plates on one node.

- **"16 configurations across three rigs" now says that all three rigs are
  Texas A&M.** #401 rejected a generalised `G_bar/G` correlation because
  69% of the variance in the pooled population falls BETWEEN rigs rather
  than within them. True, and the phrase reads as more independent than it
  is: 91-GT-3 and CR-3837 are Han's own laboratory, `lau1990` is Lau's in
  the same department. It cuts both ways and both are now recorded where
  the number lives:

  - **the conclusion is STRONGER** -- three rigs inside one laboratory
    tradition separating this far makes the offset rig-level (geometry,
    instrumentation, data reduction) rather than an artefact of comparing
    institutions, so there is even less for an angle or shape term to
    reach;
  - **the MAE 8.3% is OPTIMISTIC** -- this is not a sample of the world's
    rigs, so the 29% spread is a LOWER BOUND and Han's 1.2 may do worse
    against a non-TAMU rig.

  Nothing computed changes. Added to `runner.py`'s decision comment, the
  91-GT-3 extraction record, and the test that pins the evidence -- which
  now also fails if a fourth source joins the pool, since that is when the
  "lower bound" caveat becomes testable rather than stated.

- **The film runner reads both `x_over_sx` and `x_over_D`.** `andrei2014`
  plots streamwise distance in pitches and `murray2018` in diameters;
  treating one as the other would misplace every row by the pitch ratio,
  which looks like a modest bias rather than a bug. The row comb step and
  the row-to-sample distance now follow the declared abscissa.

- **Murray & Ireland (2018) read, and it separates #420's two error
  sources.** #420 found Baldauf + Sellers failing on `andrei2014` with the
  error changing sign (+61% at BR 1, -81% at BR 3) and could not say whether
  the superposition or the per-row closure was at fault. Murray applies
  Sellers to a single-hole **CFD** result -- no correlation involved -- and
  measures a **~2x over-prediction at M ~ 1** for a 5.75D pitch, growing with
  blowing, caused by streamwise jet interaction. Tripling the streamwise
  pitch to 17.25D restores agreement to under 10% by the third row.

  So: superposition over-predicts and grows with blowing and tightness; the
  Baldauf closure under-predicts above M ~ 1 through lift-off. `andrei2014`
  shows both, which is why its error changes sign, and its streamwise pitch
  of 9.15D sits between Murray's failing and working cases as the +61% at
  BR 1 suggests it should.

  This **qualifies** #420's "Gao's alpha cannot be fitted" conclusion. That
  holds for `andrei2014`, but the reason is narrower than recorded: alpha is
  exactly the right shape for the superposition error, and Murray's
  independent CFD confirms that error is real and one-signed. What
  disqualifies Andrei is the closure contaminating BR 2-3. The test docstring
  now says so.

  No data extracted yet. Recorded in
  `validation/cooling/extractions/murray_ireland_2018_effusion_superposition.md`,
  with a per-pitch map of what the paper actually offers: the superposition
  comparison exists **only at 5.75D** (Fig. 6) and at 17.25D streamwise
  (Fig. 9a), while the 3.0D data is contours with no superposition column.
  Figure 6 is the target because it puts measurement, CFD reference and
  superposition on one axis, and its M 0.19-0.96 sits entirely below the
  lift-off peak that wrecked `andrei2014` -- but its s/D of 5.75 is still
  outside Baldauf's 2-5, by 15% against Andrei's 47%. Also recorded: the
  measurement band to use (the paper's stated 15%) and the warning that
  every figure page is raster, so this needs pixel digitisation rather than
  the exact vector extraction `andrei2014` allowed.

- **`test_thermo_transport` no longer captures stderr**, so the C++ suite runs
  where `/tmp` is not writable. googletest's `CaptureStderr()` mkstemps into a
  hardcoded `"/tmp/"` (`gtest-port.cc`) that no environment variable
  redirects, and a failure there is fatal -- it aborted the whole 175-test
  binary, not just the two tests using it. They now read combaero's own
  warning handler through a small RAII `WarningCapture`, which is also the
  more direct assertion: that the warning was emitted, not that it reached a
  particular stream. Local `ctest` goes 35/36 to 36/36.

- **`McGreehanSchotsch1988` now warns below its `re_min = 1e4` floor**,
  naming Lichtarowicz and Idelchik as the valid alternatives. Below the floor
  the chain holds Re at the floor, so Cd stops responding to flow entirely --
  documented and deliberate, since Eq. (8) diverges there, but silently
  returning a frozen Cd is how a plenum-fed hole came to read +21.7% high.
  Warns rather than refuses, because a network can transit low Re during
  Newton iteration.

- **`EffusionPlateElement`** (#387) -- a multi-perforated wall panel as a
  network element, flow side. Geometry is given the way a plate is designed
  (pitch, hole diameter, wall thickness, inclination angle); the hole count
  follows from the panel area and is rounded, with `hole_count_exact` and
  `pitch_actual` recording what the rounding did. An inclined hole is
  correctly longer than the wall is thick (`L = t/sin(alpha)`), which is why
  effusion holes are inclined at all. Geometry reproduces Andrews et al.
  (1988) ASME 88-GT-290 Table 1: `N = 4306 m^-2` and `A/A_h` to 1.5%.

  Defaults to `IdelchikThick` -- diagram 4-18a is a thick-walled hole in a
  large wall between two plena, which is exactly a plenum-fed effusion
  plate. ISO 5167 correlations are **refused**: a panel has no pipe and no
  `beta`, so they are undefined rather than merely inaccurate. For a
  channel-fed panel the approach flow is a crossflow, so use
  `'McGreehanSchotsch'`.

  **The effusion CHANNEL needs no new element.** Because the solver balances
  mass over every element at a node, hanging a panel off each segment of a
  coolant channel gives `m_channel_in = m_channel_out + m_effusion` -- the
  coolant flow falling along the wall. That is van de Noort and Ireland's
  (2022) flow-network structure, and the resolution is the caller's choice of
  segment count.

  A single panel homogenises, and `diagnostics()` says so where it matters:
  `is_ingesting` flags a reversed drive (hot-gas ingestion, which a
  homogenised panel would otherwise average into a healthy net outflow), and
  `dP_drive` compared across panels shows whether one panel is smearing a
  gradient that deserves several. Also reports `G_coolant`, the coolant mass
  flow per unit plate area that Andrews correlates on.

  Thermal modelling is deliberately not included yet: the effectiveness split
  between internal throat convection and the external film has to be stated
  explicitly or the two double-count, and that needs its own extraction.


- **`CdCorrelation` renamed to `MeteringCdCorrelation`** and narrowed to the
  normed ISO 5167 device it actually describes:
  `ReaderHarrisGallagher`, `Stolz`, `Miller`, `Constant`, `UserFunction`.
  `Constant` and `UserFunction` are now reachable from Python, and all three
  dispatch sites (`solve_orifice_mdot`, `orifice_flow`, `make_correlation`)
  go through one factory instead of three divergent switches.


- **`Cd_McGreehanSchotsch`, `mcgreehan_schotsch::cd` and `cd_with_crossflow`
  take an `eps` argument** (#383), defaulting to `rv_smooth_eps`. Passing
  `eps = 0` recovers Eq. (17) exactly at every entry point, as it already did
  for the expansion factor. At the default the zero-crossflow value sits
  ~4e-4 above the unregularised one -- 2% of the correlation's own data
  scatter -- stated rather than silent; the paper's printed 0.60 baseline is
  still reproduced to 5e-4 at `eps = 0`.

- **Han's `G_bar = 1.2 G` is now applied at every rib angle, as published**
  (#339 plan 1b). The runner refused anything off 90 degrees on the grounds
  that 1.2 is a 90 degree result. It is -- but refusing withheld a number
  Han does publish, and turned a known accuracy limit into a missing
  answer. Faithful implementation uses the published constant; what it
  costs is measured and labelled (see docs/VALIDATION_POLICY.md).

  **A generalised replacement was considered and rejected on evidence.**
  Pooling every closed-form measurement available -- 91-GT-3 table 3
  (extracted and two-channel verified for this), Lau table 2, and
  CR-3837's per-run `Nu(R)`/`Nu(AV)` split -- gives 16 configurations
  across three rigs spanning **1.096 to 1.413** at `e+` = 300, and
  **69% of that variance is BETWEEN RIGS**:

  | source | n | mean | internal spread |
  |---|---|---|---|
  | 91-GT-3 | 4 | 1.160 | 12% |
  | CR-3837 | 5 | 1.323 | 3% |
  | lau1990 | 7 | 1.327 | 14% |

  No angle or shape term can reach a rig offset, so a "better" correlation
  would be fitting rig identity. Han's 1.2 scores MAE 8.3% against that
  population -- a defensible one-number answer. Pinned by a test that fails
  if the between-rig share drops, which is the condition under which
  revisiting would be worthwhile.

  Even at 90 degrees the measured ratio runs 1.29 down to 1.13 across
  `e+`, so the constant was always a fit, including where it was trusted
  most.


- **The cooling scorecard now aggregates all three runners** (#333).
  `validation/cooling/scorecard.py` knew only `runner.py`, so the Florschuetz
  and orifice series read as "not scored by any set" and there was no single
  place to see whether cooling as a whole was healthy. `run_dataset()` gives
  every series to exactly one runner and a rollup reports per correlation
  set -- six sets, 680 points, in one view for the first time.

  Dispatch is on an explicit `owns()` predicate per runner, not on "the runner
  returned something". Every runner deliberately returns reason-carrying
  records for series it cannot score, precisely so a silent zero never
  masquerades as a good score -- which makes truthiness useless as an
  ownership test, and was the first version's bug: the jet-array runner
  claimed the orifice series and reported them unscored.

  Rohde's series are rolled up **per velocity-head-ratio band**, never pooled:
  the total-to-static factor is 1.01 at VHR 50 and 1.33 at VHR 2, so one
  number would mix a near-exact comparison with one dominated by the
  conversion.

- **`CLAUDE.md`'s explicit-includes rule names both strictnesses.** It said
  "macOS-only implicit includes break Linux CI", which is one direction of a
  two-directional problem and misleading: a `std::max({a,b,c})` here passed on
  macOS **and** Linux and was caught only by MinGW. libstdc++ is strict about
  transitive headers; MSVC is strict about POSIX-ish macros like `M_PI`. Any
  platform can be the permissive one, and two agreeing is not evidence.

- **`ChannelElement` supports smooth surfaces only** for now. The
  user-set `Nu_multiplier` and `f_multiplier` on `ConvectiveSurface` are
  **unaffected** -- they encode nothing, default to 1.0, and remain the
  supported way to match measured data at an operating point.

- **The GUI rejects saved networks carrying removed surface types** rather than
  silently substituting smooth, which would answer a different question than
  the one asked. Ribbed reports that it returns in a later release; the other
  three report removal.

- **`ChannelElement.residuals` no longer restates pipe friction in Python.**
  It computed a `f_base` from a hardcoded Haaland/Petukhov branch, plus the
  density, viscosity, velocity and Reynolds number feeding it, on every
  residual evaluation -- including a `complete_state` call. After the pin-fin
  and dimple fixes, all of it fed nothing but a guard: `rib_friction_multiplier`
  takes no Reynolds number and `dimple_friction_multiplier` ignores the one it
  is given. Removed, along with the air-at-STP fallbacks (`rho = 1.2`,
  `mu = 1.8e-5`) it substituted whenever a state looked unphysical -- fabricated
  properties a mid-Newton iterate would silently pick up.

  Verified behaviour-neutral: residual and every Jacobian entry **identical
  across 108 cases** -- four surface models, nine mass flows including zero,
  1e-12 and reversed, three temperatures.

  Reynolds number for the dimple correlation is now computed by a small helper.
  Density cancels out of it (`rho * v` is `m_dot / area`), so no density guard
  is needed at all.

- **`ChannelElement` carries its provenance.** The docstring was one line with
  no literature. It now names the formulations, points at the files that own
  each correlation rather than transcribing them, and records the known
  Jacobian gap.

- **`test_channel_jacobians.py` checks derivatives rather than field
  presence.** Its assertions were `assert result.ddP_dmdot != 0.0`, which says
  a field was written, not that it is right -- an element-level defect 10.5%
  in size survived it. The ribbed and pin-fin cases now verify `ddP_dvelocity`
  against a central difference of `dP` in velocity.

- **A localised array's pressure drop is now the correlation's own.** The
  pin-fin and impingement correlations return `dP` directly, and
  `ChannelElement` takes it instead of converting it into a multiplier on pipe
  friction. The former route divided by a friction factor restated in Python
  and let the C++ one multiply it back in, which cancels only when both pick
  the same correlation. The Python side always used Haaland (rough) or
  Petukhov (smooth) regardless of the element's `friction_model`, so the drop
  moved by **+0.48%** for `colebrook`/`serghides` and **-24.7%** for
  `petukhov`. For the default `haaland` the drop is unchanged bit for bit.

  A channel carrying one of these surfaces now takes the array drop in both
  regimes. The array correlations are incompressible by construction, so
  `regime="compressible"` no longer layers Fanno friction on top of a drop the
  array already owns.

- **`basis_of` matched only the source folder name and its `after` field,
  so a correlation set keyed on a paper NUMBER could not recognise its own
  data.** `andrews_1986_effusion_internal` is keyed "86-GT-225" --
  deliberately, to stop it matching the 1988 study from the same lab --
  and that string appears only in the citation. The paper's own Figure 8
  therefore reported "accuracy", understating a fidelity result: the exact
  mirror of the error the key was introduced to prevent. The source
  citation now joins the haystack, and
  `test_every_scored_series_basis_is_what_its_provenance_says` pins the
  label of every scored series so the wider match cannot claim a false
  fidelity.

- **Saving and restoring the warning handler nested without bound, and
  eventually segfaulted.** `get_warning_handler()` wrapped the current
  handler in a fresh `py::cpp_function` on every call and
  `set_warning_handler()` wrapped that again, so each save-and-restore cycle
  added a layer. The layers cost nothing until a warning is emitted, at
  which point `warn()` recurses through all of them. Measured: **6000 cycles
  survived, 12000 took the process down** with no catchable Python error.

  Not exotic usage -- `suppress_warnings()` performs exactly this cycle, and
  a validation runner that probed a correlation's envelope once per sample
  reached ~6000 calls and segfaulted the whole test suite.

  The handler installed from Python is now a named `PyWarningHandler` rather
  than a lambda, so `get_warning_handler()` can recognise it through
  `std::function::target` and hand back **the original object** instead of a
  new wrapper. The cycle is idempotent: `get_warning_handler() is handler`.
  50000 cycles now pass, and the handler fires once rather than N times.

  A C++-side handler (the default included) still has to be wrapped, since
  there is no Python object to return -- but that wrapper is not a
  `PyWarningHandler`, so the chain stops at one layer instead of growing.
  Covered separately.

  Distinct from the at-exit crash fixed earlier in this release: that was a
  Python reference released after `Py_Finalize()`, this is unbounded
  recursion during the run.

- **30 functions advertised a phantom overload.** pybind11 tries overloads
  in registration order, so a second `m.def` with an identical signature can
  never be selected -- but it is still documented. `help(combaero.bernoulli_P2)`
  printed "Overloaded function" and listed two numerically identical
  signatures, inviting the reader to work out a difference that does not
  exist. 31 unreachable registrations removed across 30 names
  (`channel_mdot` had two).

  Most arose from a correlation being bound once in a general section and
  again in a domain-specific one. Nothing ever failed, which is how they
  accumulated.

  Where the two differed only in documentation, the **better-documented one
  was kept** and the other dropped, so this is a pure deletion with no
  docstring edited. That recovered rather than lost text in four cases:
  `channel_roughness` keeps its 17-material list, references and worked
  example; `channel_pressure_drop`, `fanno_max_length` and `mole_to_mass`
  likewise.

  The five genuine overload pairs are untouched -- `orifice_flow_thermo(..,
  Cd)` against `(.., cd_fn)`, `channel_flow_rough` with and without
  `k_loss_fn`, and `combustion_state`/`combustion_state_from_streams`/
  `equivalence_ratio_mass`. Guarded by
  `test_units_sync.py::test_no_phantom_overloads`, which compares signature
  strings at runtime so it cannot drift from the binding.

- **Baldauf's film-effectiveness correlation never checked its own validated
  envelope.** The bounds have been `constexpr` in `cooling_correlations.h`
  since the correlation landed -- `M_min/M_max`, `P_min/P_max`,
  `s_over_D_min/max`, `alpha_deg_min/max`, `Tu_min/Tu_max` -- under a comment
  saying that outside them the correlation is extrapolation. Nothing read
  them. `film_effectiveness_baldauf_2002` answered a hole spacing of 7.4 as
  confidently as one of 3.0, silently.

  Both entry points now take an optional `CorrelationStatus *status`,
  following `nusselt_dittus_boelter`'s contract: pass one to take the
  extrapolation flag silently, leave it null to get a warning per
  out-of-range input. **Still reported, not enforced** -- the value outside
  the envelope is bit-identical to before, because a network solve transits
  odd states during Newton iteration and refusing there would break
  convergence rather than protect anyone.

  This matters beyond a caller's own judgement: the validation harness
  derives its `extrapolated` column from this signal, so a silent
  correlation makes extrapolation indistinguishable from model error in the
  scorecard -- which is #389's fidelity/accuracy confusion in concrete form.
  Measured on the `andrei2014` rig, **all six of its (BR, DR) curves are
  outside the envelope**: `s/D` 7.37 against a maximum of 5 on every curve,
  `P` 1.0 against a minimum of 1.2 on three, `M` 3 against a maximum of 2.5
  on two.

- **A Python warning handler outlived the interpreter and segfaulted at
  exit.** `set_warning_handler()` stores its argument in a function-local
  static on the C++ side; when that argument is a Python callable the static
  owns a `py::object`, and a function-local static is destroyed by
  `__cxa_atexit` -- which runs *after* `Py_Finalize()`. Releasing a Python
  reference on a finalised interpreter is a segfault.

  Installing a handler and leaving it in place was enough to trigger it, so
  the binding's own documented example crashed:

  ```python
  combaero.set_warning_handler(lambda msg: None)  # exit 139 (SIGSEGV)
  with combaero.suppress_warnings(): pass         # exit 139
  ```

  The process did all its work and printed its output first, dying only
  during shutdown, so nothing in the suite noticed -- `pytest` never leaves a
  handler installed. `set_warning_handler(None)` was always clean, which is
  what localised it.

  Fixed with a module-teardown capsule in `_core.cpp` that resets the handler
  to the C++ default while the interpreter is still alive, leaving the static
  holding a plain function pointer. Covered by
  `python/tests/test_warning_handler_lifetime.py`, which runs real
  subprocesses and asserts on exit codes -- the only place this class of
  defect is visible.

- **Two composition warnings bypassed the warning handler**, so
  `suppress_warnings()` and `set_warning_handler()` could not reach them.
  `normalize_fractions` (all-zero input) and `convert_to_dry_fractions`
  (water-vapour-only input) wrote straight to `std::cerr`. Both now go through
  `warn()` like every other combaero warning, which makes them suppressible
  for the first time.

  Their prefix therefore changes from `Warning: ` to the handler's standard
  `[combaero] `; the message text is otherwise unchanged. Nothing in the repo
  matched on the old prefix.

- **`test_units_sync` was checking `__all__`, not the API.** The Hard Rule is
  that every exported Python symbol carries a `units_data.h` entry or an
  explicit `IGNORE_LIST` entry. The test enforced it by walking `cb.__all__`
  -- a hand-maintained list -- so any symbol added without also editing that
  list escaped the check entirely. **54 had**, and 14 of those had no units
  entry at all. The test now walks the real public surface (`dir()` union
  `__all__`), and a second test keeps `__all__` from falling behind it again.

  - units entries added for `pressure_loss`, `velocity_from_pressure_loss`,
    `zeta_from_Cd`, `Cd_from_zeta` and `molar_mass` (an alias of `mwmix`, so
    the sync check saw the alias name rather than the target).
  - `contextmanager`, `import_module`, `Generator` and `PackageNotFoundError`
    were reachable as `combaero.<name>` -- plain stdlib imports leaking into
    the package namespace. Now imported under underscored names. This is the
    same shape of leak as a dev-only import escaping into a wheel, and
    nothing in the repo referenced them.
  - the 50 genuine exports that were missing are now declared in `__all__`,
    grouped under a comment so the omission stays legible.
  - result aggregates (`OrificeResult`, `MassStream`, `WallCouplingResult`,
    `ChannelSolverResult`, `MomentumChamberResult`,
    `IncompressibleFlowSolution`), the `Registry`, the `vortex` submodule and
    the deprecated `standard_dry_air_composition` go to `IGNORE_LIST`, which
    is where the existing convention already puts result structs and modules.

- **`docs/API_CPP.md` was missing `film_superposition_sellers_and_gradient`
  and `film_superposition_corrected_and_gradient`** (from #413).

- **A validation series' stated uncertainty was its own model error** (#389).
  `uncertainty` is meant to be the band the SOURCE measured to, so `within`
  answers "does the model land inside the experiment's own precision". In 30
  series it had instead been set from combaero's measured error on that same
  series, which makes `within` ask whether the model is inside its own error
  -- a question it largely must answer yes to. The declared band tracked the
  measured RMS at **r = +0.9987** across `florschuetz1981`'s 27 series; the
  series whose `cross_check` reads "bias +31%, RMS 32%" carried a band of
  0.324.

  Repaired at source:

  - all 27 `florschuetz1981` series now carry **0.051**, from
    `sqrt(0.05^2 + 0.0097^2)` -- the primary paper's own +/-5% Nu measurement
    uncertainty at 95% confidence (`han_impingement.md` item 12) in
    root-sum-square with 0.97% digitisation precision. Deliberately NOT item
    23's 5.6% standard error, which is how well the CORRELATION fits its own
    data and so is a property of the model.
  - `han2012/fig4.47_R_vs_alpha` (0.105) and `fig4.47_G_vs_eplus` (0.088) now
    carry Han and Park's published **6%** (`han_ribbed.md` item 7), as every
    other scored series in that source already did. The 0.105 was that
    series' own RMS, in a cell whose `cross_check` said "Han's stated band is
    6%" in the same breath.
  - `han2012/fig4.193c_G_scatter` (0.069) now carries **no band**. That 0.069
    was `han_ribbed_high_re.md`'s own "Eq. (18) against the cloud, RMS 6.9%",
    and Rallabandi, Yang and Han (2009)'s own uncertainty has not been read.
    `within` renders `-`, which is the honest report; a band borrowed from a
    different paper's rig would not be a measurement either.

  **`florschuetz_1981_inline`'s `within` moves 58.7% -> 47.5%.** That figure
  appears in #333's progress comments and this is a correction to it, not a
  regression: no model behaviour changed, only the yardstick. Individual
  series now span 0.0% to 100.0%, where the self-derived band confined them
  to 33.3%-88.9% -- a series always sits near its own RMS, so it could
  neither fall wholly outside such a band nor wholly inside it. MAE, RMSE and
  bias are untouched.

  `python/tests/test_uncertainty_is_measurement.py` pins this four ways: the
  structural invariant (a band is a property of the instrument and the
  figure, so it cannot vary between scored series of one figure and
  quantity), the correlation check, and a pin per repaired source. The
  structural test is the one that binds -- but it cannot see the three
  `han2012` series, which are each alone in their group once unscored series
  are dropped, so those are held by an explicit pin. A numeric "band equals
  this series' own RMSE" rule was tried and rejected: with a constant 0.051
  over 27 series, two land on their own RMSE by coincidence.

  `lau1990` was checked and cleared. Its r = +0.97 is a false positive: three
  series, two bands (+/-5.8% Stanton, +/-10.9% friction, both Kline and
  McClintock from the paper), correlating only because friction happens to
  have both the larger band and the larger error.

- **The discharge-hole correlations were given the pipe Reynolds number**
  (regression from #409). `Idelchik1966*` and `McGreehanSchotsch` are
  functions of the HOLE Reynolds number, but `OrificeElement._effective_Cd`
  passed `Re_D`, built from the upstream channel diameter. Wrong by 1-8% for
  an orifice in a pipe, and by up to **38% for a plenum-fed hole**, where
  `D_up = 0` froze `Re_D` at a `1e5` fallback so the correlation stopped
  responding to flow entirely. Now computed as
  `4 m_dot_hole / (pi d_hole mu)` via an overridable `_hole_reynolds` hook,
  which `EffusionPlateElement` overrides to divide by the hole count.

- **Idelchik (1966) wall-orifice discharge coefficients** --
  `DischargeCdCorrelation::Idelchik1966{Sharp,Thick,Beveled,Rounded}`, from
  Section IV diagrams 4-17 and 4-18 (a hole in a large wall, `F1 = F2 = inf`,
  which is the effusion-plate geometry). `zeta` is referenced to the hole
  velocity and carries the full permanent loss, so `Cd = 1/sqrt(zeta)` is
  exact rather than a convention-dependent conversion. Valid **Re 25 to 1e6**
  -- three decades below McGreehan-Schotsch's `re_min = 1e4`, and the range
  small cooling holes actually operate in. No digitisation was required:
  every value is tabulated in the source and was read off the page at 400 dpi.
  Provenance in
  `validation/cooling/extractions/idelchik_1966_wall_orifice.md`.

  **Cross-source agreement, labelled cross-source:** against McGreehan and
  Schotsch (1988) the two are within 5% over the whole `L/d` range and
  **0.05% at the sharp-edged baseline** (0.5923 vs 0.5926) -- independent
  sources 22 years apart. That is an accuracy result, not a fidelity one.

- **A separate `DischargeCdCorrelation` selector**, with
  `DischargeHoleGeometry` / `DischargeHoleState` and analytic
  `discharge_cd_and_derivatives` returning
  `(Cd, dCd/dRe, dCd/d(U1_over_Vi))`. A hole in a wall has no pipe to form
  `beta` with and its `Cd` depends on `L/d`, `r/d` and the approach
  crossflow, none of which `OrificeGeometry`/`OrificeState` can express. One
  enum over both families forced every caller to pass a meaningless `D` and
  silently dropped the crossflow.

  Interpolation is monotone cubic (Fritsch-Carlson), not linear: linear would
  put a derivative jump at each of 14 table knots, and a natural spline would
  overshoot the flat tails and invent a `Cd` above the source's.


- **The `e+` path scored unrecognised quantities against `G`** (#339 plan
  1a). The branch ended in a catch-all `else: predicted = g`, so any
  `y_axis` that was not `G_bar` or `R_normalised` -- an absolute `R` among
  them -- would have been compared against a different physical quantity
  entirely. Nothing hit it only because every such series carried
  `scores: null`. Now explicit per quantity, with anything unpredictable
  returned unscored and with a reason.

- **`han_park_1988_angled` scored every angled series at 90 degrees**
  (#391). On the `e+` path `run_series` copied `e_D`, `p_e` and `W_H` out of
  `series.geometry` but never `alpha_deg`, a top-level field, so
  `_probe_geometry`'s default of 90 survived. Every angled series was scored
  as though its ribs were transverse -- against the one correlation set
  whose entire subject is rib angle.

  This invalidates the "RMS 8.5%" cross-check against figure 4.51's
  parallel-rib classes recorded in `han_ribbed.md`. Re-scored,
  `fig4.51_G_60par` reports MAE 6.8% and bias -6.8%, and the set moves from
  bias +4.4% to +1.5%. Nothing failed before the fix; the numbers simply
  meant something other than what they said.

- **The `fig4.46` class-label dispute is closed** (#393). The `e/D` 0.047,
  `P/e` 10, `W/H` 2 series was marked `disputed` because its `G` and `R`
  panels paired on only 2 of 4 marks and the digitised data could not say
  which symbol the marks belonged to.

  Resolved by going behind the figure rather than re-reading it. **NASA
  CR-4015** (= AVSCOM 86-C-25, Han/Park/Ibrahim, contract NAS3-24227,
  September 1986) is the report Fig. 4.46's "This study" classes come from,
  traced through Han and Park (1988)'s acknowledgement. Its appendix
  tabulates every run with a self-identifying header, so the class is found
  by label and the overplotted cluster never has to be resolved: five runs
  at `Re` 10111-64193 give `e+` 80.7, 150.5, 256.2, 485.8, 512.1, and every
  mark in both panels lands on one of them.

  Two recorded claims were wrong. The class is not "absent near `e+` ~ 80"
  -- run 22 sits at 80.7 and the `R` panel picked it up to 0.1% -- and its
  predicted ceiling of 479 was too low, the true maximum being 512.1. The
  non-pairing was overplotting, not evidence about the label.

  Both series move to `class_confidence: confirmed`, and a new test pins the
  finding rather than the prose: every digitised mark must land on one of
  the five documented runs.

- **Lau metadata quoted a withdrawn resolution** (#392). `lau1990`'s
  `metadata.yaml` and extraction doc stated that Han's `G_bar` is the
  Prandtl-normalised `G Pr^-0.57` and "not a four-wall average", citing
  `han_ribbed.md` item 10. That is item 10's **superseded** form; its current
  resolution is `G_bar = 1.2 G`, a four-wall average, with the Prandtl
  reading recorded there as withdrawn. The inverted claim was the stated
  reason for withholding Lau's `Gbar` and for a test asserting it stay out of
  the dataset.

  Confirmed from primary raw data rather than re-reading the figure: NASA
  CR-3837's appendix prints `Nu(R)`, `Nu(S)` and `Nu(AV)` per run, and
  `Nu(AV)` is their two-way mean to 0.084% over 33 rows. No code defect --
  `runner.py`'s `G_BAR_OVER_G = 1.2` was correct throughout. The test that
  pinned the wrong conclusion is replaced by one pinning the guard that
  actually matters: the `G_bar` path stays restricted to 90 deg ribs.

- **`uncertainty` now means one thing across the cooling dataset** (#333).
  `han2012` and `florschuetz1981` declare it as the band model agreement is
  judged against (3-8%); the two orifice sources added in #381 declared
  digitisation precision instead (0.2%), and for McGreehan-Schotsch it was an
  *x-axis* calibration residual being compared against a y-error. Both are now
  `null` -- neither source states an accuracy band -- with the precision
  figures kept in the metadata headers where they cannot be mistaken for
  tolerances. A series with no stated band now reports `-` for `within`
  rather than a structural `0%`.


- **Ribbed channels' pressure-drop Jacobian had no temperature or pressure
  column** (issue #378). `ChannelElement._ribbed_residuals` computed `rho`
  via `_safe_rho(state_in.density())` and never differentiated it, even
  though `dP` depends on upstream temperature and static pressure through
  `rho(T, P)` -- its own docstring already flagged this as a known,
  deliberate simplification, but it went unmeasured because
  `tests/test_channel_element_jacobian_fd.py` (built for exactly this class
  of gap) never actually added `"ribbed"` to its surface matrix, despite
  that file's own docstring earmarking it for #334.

  Adding `"ribbed"` there while extending the file for #337's impingement
  elements surfaced the gap with a number: the analytic `A.T` column was
  exactly `0.0` against a finite-difference value of roughly `-486` for a
  representative case. Not a wrong converged answer -- the residual itself
  is exact at every iterate -- but a real, non-negligible missing term in
  the Newton step direction, costing convergence speed and robustness for
  any network solving upstream temperature near a ribbed channel.

  Fixed with `density_and_jacobians`' analytic `d(rho)/dT`/`d(rho)/dP`
  (already used elsewhere in this file), chained through `_safe_rho`'s own
  floor derivative and `d(dP)/d(rho) = -dP/rho`. Composition (`Y`)
  sensitivity through molecular weight remains an accepted gap -- no
  correlation-level Jacobian for it exists yet, matching every other
  channel-friction path in this file. Verified against finite differences,
  including a dedicated static-pressure check not folded into the shared
  test matrix (that column isn't verified for the other, C++-backed surface
  types).

- **The Fanno march refines its step near the sonic point.** Fixing the
  gradient (below) gave the march a `1/(1 - M^2)` term that a fixed-step RK4
  cannot resolve: the step that is ample over most of a duct steps straight
  past the singularity in the last few percent of one that chokes. The
  reported outlet then snapped to the march grid, so the pressure drop moved
  with `n_steps` and its gradient with respect to mass flow developed a
  sawtooth -- 988% peak-to-peak at 100 steps.

  Both marches (constant-`f` and rough-wall) now halve the step whenever a
  trial step would close more than a quarter of the remaining distance to
  sonic, and relax back toward the nominal step once clear. `n_steps` sets
  the nominal step rather than a fixed one. A choked duct's outlet pressure
  is now identical to within a millipascal from 100 to 6400 nominal steps,
  where it previously moved 1.5 kPa between coarse and fine.

  Refining everywhere would have cost the same accuracy for every duct that
  never approaches sonic, which is most of them. The extra steps are spent
  only where the gradient is actually stiff. See issue #363.


- **The Fanno march was integrating an incompressible gradient.**
  `dpdx_fanno` returned `-f/(2D) rho u^2` -- Darcy-Weisbach. The Fanno static
  pressure gradient carries `(1 + (g-1)M^2)/(1 - M^2)`, and the `1/(1 - M^2)`
  term is what produces the choking singularity. Without it the flow never
  choked however long the duct: marched over its OWN `fanno_max_length`, where
  M must reach 1 by definition, it reached M = 0.37 / 0.55 / 0.72 from inlet
  Mach 0.3 / 0.5 / 0.7 and reported `choked = false` every time. It now agrees
  with an independent integration of the standard ODE (M = 0.9841 against
  0.9714 at L*).

  `fanno_max_length` returned **half** the correct length: its estimate divided
  by `4f` where `f` is Darcy, putting it at L*/4, and with a bracket of only
  twice the estimate the bisection saturated at its own upper bound. It now
  matches the tables to 0.02% from M = 0.3 to 0.85.

- **The compressible channel confused static and total pressure, three ways.**
  It evaluated density at `P_total` ("uses P_total as proxy for P_static
  (low-Mach approximation)"), so at M ~ 0.76 it was 1.59x too dense and
  reported converged solutions for flows 1.37x past its own choked limit. The
  inlet static state is now inverted from `(Pt, Tt, m_dot, A)` -- energy,
  isentropic static/total, continuity -- as a 1-D root find on Mach, using the
  entropy-based relations `stagnation.h` already provided.

  The march returned a STATIC drop while `ChannelElement` consumed it as a
  stagnation drop; in Fanno flow the flow accelerates and converts static head
  into dynamic, so the two differ by 1.09 - 1.58x over M = 0.34 - 0.64. It now
  returns `Pt_in - Pt_out`.

  The Jacobian chained through the march's inputs as if they were the element's
  unknowns, which held only while density came from the total state. Chained
  through the inversion instead, the error against a central difference fell
  from 4.7e-2 to 3.1e-8.

### Deprecated

- **`VortexElement(r_c=..., n=...)` (#661).** The free vortex has no Vatistas core, so both arguments are ignored with a `DeprecationWarning`; they are removed in 0.8.0. The `vatistas_*` functions and `VatistasVortex` are unchanged.

- **`BorderCarnotLossElement`, for removal in 0.8.0.** Its residual is exactly a constant head loss.
  - **Migration:** `PressureLossElement(id, from_node, to_node, correlation=ConstantHeadLoss(zeta=4 * (1 - cos(0.75 * delta)) ** 2, area=area), area=area)` reproduces it to round-off (1e-11 measured), with `delta` in radians.
  - **Why it goes:** its only addition is that zeta formula. It was derived for a junction's lateral branch, and `MultiPortChamberElement` and `ConstantKTeeElement` already carry that loss, so adding it to their ports double-counts (#272). It was never compared against bend data.
  - **Usage:** nothing in the package, GUI or examples uses it. Constructing one now emits a `DeprecationWarning`.
  - **C++:** the backing `border_carnot_loss_residual_and_jacobian` and `border_carnot_L` go with it.

- **Every deprecated name now warns and names its removal release, 0.8.0.** Before this, four were deprecated only in a comment or a docstring, so a user had no notice before the name disappeared:
  - `standard_dry_air_composition()`: use `species.dry_air()`.
  - `smr_wgs_equilibrium` / `smr_wgs_equilibrium_adiabatic`: use `reforming_equilibrium` / `reforming_equilibrium_adiabatic`. These are equivalent replacements only since the simultaneous-equilibrium fix; before it, the replacement was 22 K off at a 900 K pre-reformer inlet.
  - `combaero.network.ConstantLoss` / `LinearThetaLoss`: use `ConstantFractionLoss` / `LinearThetaFractionLoss`. They warn on access and are no longer star-exported.
  - `init_strategy='incompressible_warmstart'` already warned and now names 0.8.0 too.
- `python/tests/test_deprecations.py` is the inventory: each deprecated name must warn, name 0.8.0, and agree with its replacement.

### Removed

- **The solver-derived junction barrier weight (#272).** Removed: `scaled_penalty_alpha`, `BARRIER_SLACK_FRACTION`, `MultiPortChamberElement.effective_penalty_alpha()` and the solver hand-off.
  - **Why:** with the closure-consistent seed the soft barrier is seldom reached. Derived and fixed weights now give identical outcomes over 600 scaled random junctions (sizes 1e-6 to 1e4) and the full junction scorecard. The control test's documented exit was deletion.
  - **What remains:** `soft_penalty_alpha` (default 1e11) is used directly. Set it explicitly for a network far from ~0.1 kg/s that parks a port in the barrier.

- **The area change's Mach correction k1 = 1 + 0.35 M^2 (#731).** Removed: the `Mach` argument of `sharp_area_change` and `conical_area_change`, `MACH_CLAMP`, `MACH_ACCURACY_LIMIT`, `AreaChangeResult.mach_clamped`, and the solver wrappers `area_change_residuals_and_jacobian`, `conical_area_change_residuals_and_jacobian` with `AreaChangeElementResult`.
  - **Why:** k1 had no source, and the network never applied it (its residuals passed Mach 0). Compressibility is now in the closure itself (see Changed).
  - **Migration:** drop `Mach=`. For the network form use `area_change_flow` / `area_change_element_flow`. The diagnostic `dP_static` is now `dPt_loss`: it was always a stagnation-pressure loss.

- **The 'Miller' orifice correlation.** No published source matched the C++ formula, `0.596 + 0.031 beta^2 + 0.5 beta^2.5 / sqrt(Re_D)`. The network's 'Miller' had been computing Stolz.
  - `MeteringCdCorrelation.Miller`, `orifice::Cd_Miller`, the GUI option and the schema value are gone.
  - A network or saved GUI file that names 'Miller' is refused with a message: pick 'Stolz' or 'ReaderHarrisGallagher'.

- **`channel_compressible_residuals_and_jacobian`, `channel_compressible_mdot_and_jacobian`** and the barrier constants `kChannelChokeMinMarched`, `kChannelChokeBarrierKappa`, `kChannelChokeBarrierBeta` (#481). They posed the channel as "given m, find the drop", which has no answer past choke; replaced by `fanno_channel_flow`. Their unit tests went with them; the march's own grid-independence test stays.

- **`CombustorNode.set_fuel_boundary`** (#481). It stored a boundary nothing read: fuel only ever counted when wired through an element. Use an element from a MassFlowBoundary, or an injecting MassFlowBoundary.

- **`Cd_thick_plate`, `Cd_rounded_entry`, `Cd_orifice`,
  `orifice::thickness_correction`, `orifice::Cd_rounded`, and the
  `IdelchikThick` / `IdelchikRounded` / `BohlThick` / `BohlRounded` members
  of the Cd correlation enum.** Replaced by the Idelchik members above.

  These computed an ISO 5167 `Cd` and multiplied it by a correction factor.
  The base does not apply -- a thick-edged or rounded orifice is not the
  normed device -- and the factor was an unlabelled fit to Idelchik diagram
  4-12b that had lost its tail: within 5% to `r/Dh = 0.06`, then **-27% at
  0.08 and -78% at 0.12**. The thick-plate factor carried a `5.67`
  commented "friction loss calibration factor" and a `[0.5, 1.3]` clamp.

  Two defects went with them:
  - an **11.1% jump discontinuity in `Cd` at `Re_D = 1e5`** (0.881391 ->
    0.979326), a hard C0 break in a solver input;
  - `Cd_rounded_entry` at `r/d = 0` **silently returned `Cd_Stolz`**.

  None of the 16 tests over that code would have caught either: all were
  directional. `BohlThick`/`BohlRounded` fell through to the Idelchik
  implementations, silently returning a different correlation than the one
  asked for; Bohl was dropped rather than implemented because its orifice
  content (DIN 1952 normed orifices, standard loss coefficients) is already
  held by ISO 5167 and Idelchik, so it plugs no coverage gap.

- **`OrificeElement`'s `'Auto'` correlation arm**, and the `'ThickPlate'` /
  `'RoundedEntry'` strings with their GUI dropdown entries. `'Auto'` picked a
  correlation from `r` and `t` behind the caller's back, which is how a
  rounded-entry request came back as Stolz. The new strings are
  `'IdelchikThick'`, `'IdelchikBeveled'`, `'IdelchikRounded'` and
  `'McGreehanSchotsch'`; an unknown correlation now raises rather than
  falling through to auto-selection.

- **Analytic `(f, J)` for the McGreehan and Schotsch discharge coefficient**
  (#383). `mcgreehan_schotsch_1988_cd_and_derivatives` returns
  `(Cd, dCd/dRe, dCd/d(U1/Vi))` via forward-mode `DualN`, cross-checked
  against finite differences to 3.7e-7 relative. `r/d` and `L/d` are geometry
  and carry no partials. This is what `OrificeElement` needs before the
  correlation can be selected there -- without it the Jacobian would be
  knowingly incomplete, since `Cd` moves 0.60 to 0.99 over the geometry range
  and is halved by crossflow.

  **The two hazards needed different treatments, and measurement decided
  which.** Eq. (17)'s `Rv^0.6` and `Rv^0.9` diverge in slope at `U1/Vi = 0` --
  the default, and the physically right value for a plenum-fed jet plate -- so
  the input is **regularised in the model**: `U1/Vi -> sqrt((U1/Vi)^2 + eps^2)`
  with `rv_smooth_eps = 3e-5`. Below `Re = 1e4` the value is held at the
  validity floor and the derivative was exactly zero; softening the *value*
  there was measured and **rejected** (a soft-max recovers only 16% of the
  live slope for 7.1e-3 on `Cd`, in a region the correlation does not cover),
  so the **Jacobian alone** is continued from the floor. No reported `Cd`
  changes for that one, and `dCd/dRe` becomes continuous across the floor.

  `eps = 3e-5` sits between two bisected walls, both sourced: worst
  `|dCd/du|` within 10x the physical derivative scale (0.663, measured on the
  exact chain over the range Figs. 4-6 carry data for) gives `eps >= 1.86e-5`;
  departure from Eq. (17) within 10% of the `+/-0.02` scatter the correlation
  sits in (**the paper states no error statistic at all**) gives
  `eps <= 4.15e-4`. The default sits 15% across, deliberately towards the
  low-smoothing end -- over-smoothing stalls this solver harder than a stiff
  Jacobian does.

- **Every cooling correlation whose provenance did not survive review.** A
  correlation-by-correlation audit against the cited sources found that the
  base convective layer is exact and almost nothing above it is. Removed:
  `rib_enhancement_factor`, `rib_friction_multiplier`, their `_high_re`
  variants, `dimple_nusselt_enhancement`, `dimple_friction_multiplier`,
  `pin_fin_nusselt`, `pin_fin_friction`, `impingement_nusselt`,
  `film_cooling_effectiveness` and its `_avg` and Sellers multi-row
  companions, `effusion_effectiveness`, `effusion_discharge_coefficient`,
  the `channel_ribbed`/`channel_dimpled`/`channel_pin_fin`/`channel_impingement`
  wrappers, the seven finite-difference `_and_jacobian` helpers built on them,
  their pybind11 bindings and `units_data.h` entries, the `RibbedModel`,
  `DimpledModel`, `PinFinModel` and `ImpingementModel` surface types, and the
  matching GUI node types.

  Why removal rather than repair: repair needs a known-good target and the
  citations did not provide one. `rib_friction_multiplier` returned **1.4534**
  where the only rib datum in the repository gives ~6.3 -- a factor of 4-5, in
  the pair the element actually shipped. The dimple pair swept depth and
  spacing that its cited source held fixed, and its `S_d` validation would have
  rejected that source's own geometry. `impingement_nusselt` cited a
  correlation whose defining term is the crossflow-to-jet mass flux ratio and
  had no mass flux argument. `film_cooling_effectiveness` decayed
  exponentially where film effectiveness decays as a power law.
  `effusion_effectiveness` used `I = M^2*DR` for a momentum flux ratio that is
  `M^2/DR` -- with the correct definition present but unused in the same file.

  The tests did not catch any of it because they measured the code against
  itself: `assert 1.0 < multiplier < 10.0` for a value 4-5x off, and
  `assert 1.3 <= f <= 2.2` against the function's own clamp.

  Migration: pin `combaero~=0.6` to keep the previous behaviour. Rebuilt
  correlations land per issue #339, ribs first (#334).

- **`dimple_friction_multiplier_and_jacobian`** (C++, its pybind11 binding and
  the `combaero._solver_tools` re-export). It reported a derivative w.r.t.
  `Re_Dh` that was **identically zero at every Reynolds number**, because
  `dimple_friction_multiplier` accepts `Re_Dh` and never uses it -- the
  multiplier is a function of dimple geometry alone. It obtained that zero by
  central finite difference, which the Solver (f, J) rule in `CLAUDE.md`
  forbids. Nothing called it: `channel_dimpled` uses the plain multiplier, and
  the archived note claiming otherwise was wrong. Callers wanting the value use
  `dimple_friction_multiplier`; there is no derivative to want.

### Fixed

- A thermal wall's film is read against the gas before the holes' coolant mixed in, as Baldauf 2002 and Colban 2011 define the effectiveness: when the covered side is fed by the plenum or momentum chamber the holes discharge into, its T_aw shifts from the mixed gas to the node's other inflows mixed without the holes' stream, and the density ratio takes that temperature. The film read the wall cooler by about (1 - eta) times the mixed-in coolant's share of T_inf - T_c (12.6 K on a 3 %-coolant test case). New `{wall}.film_gas_undiluted` says which gas was used. Where eta is below that share (an effusion panel at high blowing ratio), a filmed wall now reads hotter than the same wall unfilmed, which assumes the coolant fully mixed (#694).

- **`VortexElement` is a free vortex (#661).** It used to add the static radial-equilibrium rise of a Vatistas vortex to the *total* pressure at constant Tt. That is neither a free vortex (no Pt change) nor rotor work (Tt rises with Pt), and it destroyed entropy for every omega > 0.
  - **Now:** Pt and Tt are conserved. r * V_theta is carried from `r_in` to `r_out`, with the inlet swirl `v_theta_in` defaulting to solid body, omega * r_in.
  - **Reported:** the static pressure at each radius, isentropic from the stagnation state (`P_swirl_in`, `P_swirl_out`, `dP_swirl_static`).
  - **Results change:** a network with a spinning vortex no longer gains total pressure across it. The geometry arguments are keyword-only, so an old positional call fails instead of being misread.
  - **GUI:** files saved before the change load as a free vortex, with a one-time note on the element. New elements default to `r_in` = 0.05 m.

- **The node type between a pipe and an orifice no longer moves the flow (#651).** With a velocity of approach the orifice drives from the upstream static state, and it read that state from its node. A `PlenumNode` has no velocity by definition, so between pipes it reported the pipe's stagnation pressure as p1, and the orifice's p1 and rho1 came out a dynamic head too high: 1.5% more flow at beta 0.7 (pipe Mach 0.18) than behind a `MomentumChamberNode`.
  - **Now:** the orifice takes the approach PIPE's static state, the node's (Pt, Tt) at the pipe's mass flux (`approach_static_pressure` in C++), so beta, E and p1 all refer to the same pipe. Plenum and momentum-chamber networks give the same flow to 1e-6.
  - **The incompressible channel** read its friction density at its node's static pressure, with the same error out of a plenum (3% in friction at Mach 0.18). It now uses its own inlet static state too. Out of a pressure boundary or plenum this lowers the density by O(Ma^2): 0.1% at Mach 0.04, 8% in the drop of the Mach 0.6 test duct, which the incompressible regime does not cover anyway.
  - **Solver:** with the approach reading the flow, an orifice with E now solves its own consistency, m* = m_calc(m*, x), as the discharge-hole Cd already did (#582). At the iterate's flow the forward kernel ran backwards wherever 0 < Pt - p2 < q_pipe, and a beta 0.89 orifice behind a plenum stalled.

- **A critical-mode ejector discharging straight into a pressure boundary now converges (#637).**
  - **What failed:** the GUI joined the outlet port's chamber to the boundary with a total-pressure-equal link. That made the reservoir recover the jet's whole head, which has no root once the chamber is sonic, so below its P_c* the ejector stalled (|F| 2.2e3 at 1.013 bar).
  - **The new link, `PortDischargeElement`,** takes the boundary at the chamber's static pressure: the head is lost, as #360 does for a channel exit. When choked, it holds the chamber at the sonic static pressure of its flux (under-expanded).
  - **Result:** the entrainment ratio now decouples below P_c*.
  - **Subcritical behaviour changes too:** a GUI ejector or tee port wired straight to a boundary loses its exit head into the boundary, so the ejector entrains less, by about 7% at 1.5 bar in the #637 case. Give the boundary `coupling="total"` to keep the old link.

- **A duct's heat transfer is evaluated at its own static state (#636).**
  `ChannelElement` (and `CrossflowSegmentElement`) and
  `PressureLossElement` passed the feed node's temperature to the surface
  correlation as the duct's static temperature. At a plenum or a chamber that
  temperature is the stagnation T, so T_aw = T + r v^2/2cp came out ABOVE the
  stagnation temperature, by about (1 + r)(Tt - T): +2.6 K at Mach 0.1,
  +22 K at 0.3, +51 K at 0.5 for a 1600 K stream (+30 K in a GUI case).
  - The static state now comes from the feed's (Tt, Pt) and the mass flux by
    the isentropic acceleration: new C++ `duct_static_state` (bound as
    `cb.duct_static_state`), the subsonic root of chamber_static_state at
    the static P whose stagnation pressure is Pt, or the sonic state above the
    choked flux. T_aw now lies between the static and stagnation temperature.
  - Past choke (a solver iterate beyond the duct's capacity) the velocity is the
    sonic one too, and h and T_aw do not move with the flow.
  - The Jacobian stays exact: h and T_aw partials are chained onto the feed's
    T and Pt and the element's flow; the wall relay gains the Pt columns.
  - New `thermo::P_FLOOR` (1 Pa): an unphysical iterate (Pt <= 0) is
    evaluated at the floor instead of returning NaN.
- **The solver's stall detector also catches a crawl near the root.** Below
  its far-from-root gate it used to stop LM only on a dead plateau. A cold
  attempt improving by less than 5% over three 5 s windows now hands over to
  the warm-start retry as well (the 5-row hot impingement array crept for
  200 s at |F| ~ 0.075; the retry solves it in 0.3 s).

- **An ejector root with a port flowing backwards is no longer reported as a solution (#619).**
  - **Why:** every regime the element models has its primary and suction streams entering the mixing plane and the mixed stream leaving the outlet. Past a jet pump's shut-off head the suction reverses, and the residual rows' continuation there has roots no ejector produces. One of them "converged" with 0.47 kg/s entering the 1.7 bar outlet, rising to a 4 bar mixing plane and leaving through the 1 bar suction.
  - **Now:** such a root is demoted to INCONSISTENT. The tolerance is 1e-6 of the largest port flow, so operation exactly at shut-off is kept.
- Stress harness: incompressible ejector draws build their port chambers in the draw's regime (they defaulted to compressible). Their outlet pressure is now drawn below the jet pump's shut-off head, 2 (Pt_p - Pt_s) b (1 - b/2) with b = A_e/A_3.

- **Junction validation harness (#272).**
  - **`random_robustness.has_root` judges a pressure-driven root by the direction of (K_straight, K_branch)**, not by the ratio K_s/K_b.
    - **Why:** the ratio's pole where K_b crosses zero hid real roots (four solvable draws in 240 called rootless) and invented one through the pole.
    - **Result:** all 7 changed verdicts in the first 400 draws now agree with the solver.
  - **imposed_q endpoints (q = 0 or 1) are scored at their one-sided limit, 1e-3 inside the range.**
    - **Before:** all 173 such records went unscored. Exactly at the endpoint a port carries no flow, and its residual kinks there.
    - **Now:** 165 are scored, and each record is flagged `endpoint_limit`.
    - **Wang:** 188 of 200 points scored, MAE 0.131 (interior 0.108).

- **Junction seed lands on the closure's operating point (#272).** `NetworkSolver` now starts a three-port `MultiPortChamberElement` from the split its own Mynard closure predicts, not a Bernoulli share of guessed pressures.
  - **Pressure guess:** a `LosslessConnectionElement` carries no pressure drop. Before, it took a share like any lossy element, which put a junction's straight port above its own supply.
  - **Pressure-driven networks** (every port at a pressure boundary): split and flow level come from the closure, taking the root that carries the most flow.
  - **Flow-driven networks** (a `MassFlowBoundary` on the common port): the imposed flow is kept, and the split is the closure root nearest the Bernoulli share.
  - **Measured on the junction scorecard (2546 records):** +38 converged, -1. The loss is an Idelchik q = 1 endpoint in a documented worst cell. Bassett three_pb now converges at lateral fractions 0.5-0.8, where it failed before.
  - **Random harness:** unchanged at 157/240.
  - **Below q = 0.5 three_pb has two genuine roots,** e.g. (0.2, 0.33). The seed takes the higher-flow one.

- **GUI solver errors name what is on the canvas (#321).**
  - A tee or ejector port is shown as the node's label and the letter it draws: `'__tee_jct__node_..._branch'` becomes "T1 (Tee) port branch (b)".
  - Every node id in the message is followed by its label and kind.
  - The two common junction errors carry a hint:
    - a port wired against the declared direction points at the merge / branch setting;
    - a port wired straight to a boundary asks for an element in between.
  - Applies to the error of a failed solve and to the solver message.
- **"mpce tee" no longer appears in user-facing text.** The probe's target list and the probe node now say "Tee". The saved type string `mpce_tee` is unchanged.
- **`scripts/run-gui.sh` uses pnpm**, as the repo does, instead of npm.

- **Every example under `python/examples/` runs again; none are skipped.**
  - `combustor.json` declares its merge chamber: the bypass orifice is a side stream, and the chamber inherits its main inlet's area. This fixes `network_runner_bc_swap.py` and `network_runner_combustion_sweep.py`. The sweep example used to exit cleanly with every row failing; it now prints its convergence count and fails if no point converges.
  - `combustor_pressure_loss_example.py` is rewritten around a `PressureLossElement` with `LinearThetaFractionLoss` on the combustor's outlet edge, the place combustor loss lives since it left `CombustorNode`. It checks the solved loss against the correlation.

- Junction Jacobians now carry each port's composition (and, for `MultiPortChamberElement`, temperature) dependence. The rows depended on it without saying so, and on a CO2 / air merge they were off by up to 7.6e-5 (scaled): `MomentumChamberNode`'s Pt closure 7.6e-5, `ConstantKTeeElement`'s common-port head 7.5e-5, `MultiPortChamberElement`'s rows through the port density 3.9e-6. All rows now match central differences to 2e-7. `momentum_chamber_residual_and_jacobian` returns `d_res_dY`.

- Junction port nodes now take their flow by its sign. A port running against its declared direction used to have the junction's outflow counted as a negative inflow, cancelling the real one, so the port node saw no flow and its Pt fell to P (0 Pa against a 711 Pa head). This affects every `MultiPortChamberBase` junction, but only at states with a reversed port: fully forward states take the unchanged path. On 400 random-BC junctions every converged Mynard-tee result is bit-identical, and three failing draws changed failure class.
- `IdealMomentumTeeElement` accepts reversed converged ports and reports d(row)/d(port T, Y), through the port density. Without them its rows were off by up to 8.7e-4 (scaled) on non-isothermal junctions.

- GUI: the tee result `mass_flow_ratio` was labelled "ṁ str/ṁ com", but every tee closure reports the branch share. It is now labelled "ṁ br/ṁ com".

- Ideal momentum junction loss factor (#526): the factor now fades to 1 as a port's flow share drops below about 1% of the throughput. With `f != 1` a receiver carrying vanishing flow sat `(f - 1)` supplied heads away from a supplier's static pressure: 4 heads at `f = 5`, across a flow change of only 1e-4 of the throughput. The solve stalled on that step, or landed in it and was refused as an "unphysical artifact root", whenever a branch shut off or reversed. A port carrying 5% of the throughput keeps 96% of the factor's effect, and 10% keeps 99%; `f = 1` is unchanged.

- **A rejected junction root no longer tells the user to try `init_strategy='analytical_pt_prop'` when that start produced it (#536).**
  - **Before:** re-solving from that start reproduced the same root.
  - **Now:** the message points at the junction's geometry instead; a port area that differs from the connected element's is the usual cause.
  - **Still suggested:** from any other start, or from a user `x0`.

- **GUI: a tee port whose channel has a different bore meets it through a sharp step (#536).**
  - **Before:** the straight port always took the main area, whatever channel was connected, with nothing joining the two.
    - With a narrower straight channel, the ideal momentum tee raised that channel's total pressure and the solve was rejected.
    - Mynard contracted losslessly.
  - **Now:** where a directly connected channel's area differs from its port's by more than 5%, the GUI inserts a hidden sharp `AreaChangeElement` (Weisbach contraction / Borda-Carnot expansion).
    - The tee's result reports each step's zeta, total-pressure drop and area ratio (`step_<port>_*`), shown under "Port steps".
    - A port with a step is modelled, so `port_area_mismatch` (#539) does not count it. The diagnostic still flags unmodelled steps in networks built in a script.
    - Solver errors name a step as `the area step at <tee> port <p>`, not by its internal id.
  - **On the reported network:** the ideal tee now solves. The ideal and Mynard tees put 0.0518 and 0.0495 kg/s down the straight leg.
- **A failed solve no longer aborts when an element cannot evaluate the diverged end state** (e.g. a negative pressure). That element reports no diagnostics, and the auto-retry runs. A converged state that cannot be evaluated still raises.

- **`AreaChangeElement` couples stagnation to stagnation pressure (#541).**
  - **Before:** it set the downstream static pressure to `Pt_up - loss`. A `MomentumChamberNode` downstream (Pt = P + q) then gained its dynamic head across the step, so Pt rose across a loss element. This node is the one the GUI inserts between two elements.
    - Example: a 0.1 -> 0.05 m sharp step at 0.06 kg/s raised Pt by about 250 Pa.
  - **Now:** Pt drops by exactly the element's own loss, as with `ChannelElement` (#175).
  - **Unchanged:** networks with a plenum or boundary downstream.

- **A discrete pressure loss follows the flow (#550).**
  - **Before:** `PressureLossElement` applied `Pt_to = (1 - xi) Pt_from` in its declared direction even when the flow reversed.
    - A reversed constant fraction raised total pressure by `1/(1 - xi)`, in solves reported successful. The second-law check found it (#545).
    - A reversed head loss lost nothing, because its correlation returns 0 for `m_dot <= 0`.
  - **Now:** the loss runs from the flow's upstream node at `|m_dot|`, blended through zero flow as the area change is.
  - **New: an optional design flow `m_ref`** for `ConstantFractionLoss` and `LinearThetaFractionLoss`. At and above it the loss is exactly the fraction; below it the loss scales with the dynamic head and vanishes at zero flow. Without it, a reversing flow must jump through the fraction's pressure difference at zero flow, which the solver cannot do. The GUI's Discrete Loss gains a "Design flow" field.

- **A converged root whose network balance does not close is reported INCONSISTENT** (#552). The `__ledger__` relative mass or energy closure must stay within `combaero.network.second_law.LEDGER_CLOSURE_TOL` = 1e-6; the suite's accepted solves close to 4.2e-8 and 1.1e-7.
- **Withdrawing more heat than a stream holds no longer aborts the solve** (#552). The mixer raised from the T(h) inverse at the starting point; it now continues below the solver floor, and the node is reported on the temperature clamp.
- **A compressible orifice fed at Tt <= 200 K flows** (#552). The isentropic solvers held the static T at 200 K, the old end of the thermo data, which zeroed the flow and gave a spurious zero-flow root. They now go down to `thermo::T_FLOOR`.

- An orifice or incompressible channel declared against its flow no longer reads the wrong side (#554). It used the declared inlet's density and drive, and the orifice also used the declared upstream pipe's velocity of approach. The same network declared the two ways disagreed: an orifice passed 0.76x the flow, and a compressible orifice declared against its flow did not converge at all. Both now evaluate the node the flow comes from, and the declaration is only a sign convention: identical results in both regimes.
- The cold-start pressure guess of an interior node is clamped to the pressure boundaries' range on networks driven by them alone. Walking the declared directions, a reversed element had seeded a node below every boundary.
- The compressible GUI default effusion liner now solves, matching the incompressible solve to 4e-6 (it was a strict xfail under #585).

- **A `MomentumChamberNode` closes on its static state (#576).** The node's `T` is the mixed stagnation temperature.
  - **Before:** the closure used that T as the static T for density, sound speed and Mach. At fixed Pt this capped a chamber at 92 % of the choked mass flux (peak at M 0.913), and Newton stalled on any chamber asked to carry more, such as an ejector outlet pipe of about the mixing-duct area (#567).
  - **Now:** the static T is solved from `h(T0) = h(T) + v^2/2`, and the closure passes up to the true choked flux. Its Jacobian is closed form.
  - **Junction ports (tee, Mynard junction, ideal tee, constant-K tee) use the same static density**, so a port's momentum rows and its node's Pt describe the same flow. Fixing only the chamber raised the Wang 2014 cross-source junction error from 0.131 to 0.194 and made it grow with Mach. Fixing both lowers it to 0.127, flat in Mach.
  - **Reporting:** a chamber's `mach`, `rho`, `a`, `velocity` and `Re` come from the static state, and `T_static` is new.
  - Remaining consumers that still read the chamber's T as static are tracked in #578.

- **Merge-chamber face, chamber heat transfer and effusion gas side read a chamber's static state (#578).** Follows #576: a `MomentumChamberNode` carries its stagnation T.
  - **Merge face:** `chamber_merge_face_state` took the chamber's T and the main stream's T as static. It now solves the face's static state from the impulse and energy equations together, and the outlet's from the chamber's own static-state solve. Its sensitivities are analytic; the temperature columns were central differences before. On a compressible liner (outlet M > 0.3), the old face missed the exact momentum balance by 1.95 N.
  - **Chamber heat transfer:** the correlation sees the static T, and the velocity from the static density.
  - **Effusion plate:** its gas-side density on a chamber is the chamber's static density.
  - **C++:** `chamber_static_state` is now public in `stagnation.h`.

- **Effusion plates and discharge-hole orifices no longer trap the solver at low flow (#582).**
  - **Cause:** the hole Cd rises steeply with the hole Reynolds number (IdelchikThick: 0.37 at Re 25, 0.86 at Re 1000). Evaluated at the solver's current flow, the element's row was non-monotone, and |F| had a false minimum at the table's Re floor. The GUI-default effusion liner stalled there in the compressible regime.
  - **Now:** the element solves its own Cd–Re consistency, m* = Cd(Re(m*)) K, and its row is m - m*. Same root and same physics; the row is linear in the element's flow.
  - **Solver:** a residual at machine precision (|F| <= 1e-12) is accepted as converged even when MINPACK's progress test misfires at the exact root.
  - **Solver:** the Levenberg-Marquardt fallback on merge-chamber networks starts from hybr's best iterate instead of the cold start.

- **Compressible crossflow and impingement chains solve their momentum balance (#585).**
  - A compressible `CrossflowSegmentElement` (and so `ImpingementCrossflowElement` and the `EffusionLiner` backside) subtracted its station momentum drop, a pressure, from the channel's mass-flow residual. The row mixed Pa with kg/s.
  - In a reservoir-bleed-outlet network at 1 kPa drive, the bleed station reversed and both segments passed about 1/50 of the incompressible flow, with a 1 um segment holding 1 kPa.
  - The station drop now shifts the Fanno discharge target, so Pt_in drives against Pt_out + dP_momentum. The two regimes now agree to 4e-4 at Ma 0.12, and the gap falls as Ma^2.
  - The default regime is incompressible, so only networks that select `regime="compressible"` change.

- **A warm-started compressible channel no longer returns a wrong flow (#585).**
  - The C++ Fanno solve (`fanno_channel_flow`) brackets the exit Mach around the caller's guess, the channel's previous exit Mach. The bracket widened only to about 3.6x the guess. A root above that was never enclosed, and the solve returned the end of the bracket.
  - A short duct warm-started at M 0.01 then passed G = 13.5 kg/(m^2 s) for every drive from 5 Pa to 400 Pa; the true flux rises from 33 to 219.
  - A guess the widening cannot bracket now falls back to the cold bracket down from M = 1. The flow no longer depends on the warm start.

- **A failed compressible solve retries from the incompressible solution whatever the network (#585).**
  - The outlet-referenced incompressible warm start rescued only networks with a junction (`MultiPortChamberBase`). A duct-fed `EffusionLiner` has none: its cold compressible solve stalled, and the same network seeded from its incompressible root converges in one pass.
  - Networks without a junction keep their whole timeout for the cold attempt, and retry only in the time left. A junction network's 40 / 25 / 35 budget split is unchanged.
  - The cold attempt's "did not converge" warning is no longer emitted when the retry runs, so a converged result no longer warns. The `RuntimeWarning` naming the retry stays.

- An incompressible channel discharging into a pressure boundary now loses its exit dynamic head, as the compressible channel always did (#360). Before, the two regimes disagreed by O(1) at any Mach: an incompressible duct's drop was about a third of the compressible one at Ma 0.025, and a jet pump behind an incompressible outlet pipe read an entrainment ratio 2.5x too high. Ribbed and pin-fin channels follow; crossflow segments, whose station nodes already hold static pressure, opt out (in both regimes, which also removes a double-counted exit head in the compressible path) (#588).
- The second-law check closes a merge chamber's main-inlet duct at the chamber's main face, which is what the duct delivers, not at the chamber's mixed total pressure, which the side jets can lift above the duct's own inlet (#588).

- A failed cold solve of a junction network that is incompressible throughout is now retried from a compressible warm start: the same network solved with its incompressible elements compressible, which shares the root at low Mach. Before, such networks got no retry at all. Reported as `compressible_warmstart` (#588).

- The ejector's subsonic jet-pump regime no longer hands its nozzle losses back. Its flow rows carried Huang's nozzle efficiencies (0.95 / 0.85), but its mixing took each stream's loss-free velocity. With suction and outlet at the same pressure, it held an entrainment ratio of 6.62 with no primary drive and zero entropy generation. Each stream now mixes at the velocity its nozzle delivers. The same case needs ~580 Pa of drive (omega 4.71), and the result agrees with the incompressible jet pump to 0.15%. The critical (double-choked) regime is unchanged (#588).

- The ideal momentum tee's energy check measures the balance volumetrically, `sum m (Pt - p_jn) / rho`, as the closure's guarantee is stated, with the port densities its rows used. Scored on mass-weighted Pt, a lossless split whose ports differ in density gained 1.1e-7 and was refused as an artifact root; it dissipates and now passes (#588).

- Species thermo beyond its NASA-9 data (below 200 K; above 6000 K for H2O, NH3, CH4 and the alkanes, above 20 000 K for the rest) no longer runs the end polynomial on into negative or runaway cp (cp/R -25 for CH4 at 50 K, -241 for H2O at 20 000 K). Under both thermo models, cp now fades C2 to a constant over 0.25 T_lo below and 0.1 T_hi above, with h and s its exact integrals, so cp stays positive and h, s and u monotone at every T. Inside the data nothing changes (#606).
- `calc_T_from_h`, `_s`, `_u`, their `_mass` forms and `calc_T_from_sv_mass` use a safeguarded Newton (bracketed, bisection fallback) and raise `RuntimeError` for a target no positive T reaches, instead of silently returning the last iterate (#606).

- **Pressure-loss, area-change, Border-Carnot and vortex rows are scaled as pressure rows** (#614). Their residuals are total-pressure balances in Pa, but the solver's row scaling filed them as mass-flow rows, which weighted them 1e5-1e9 too heavily in the norm hybr and LM minimise. On hard networks the LM fallback then walked pressures to -4e14 Pa (where the Jacobian went non-finite). A 160-network stress sweep converges 42 random networks instead of 33 (13 gained, 4 lost); in the test suite these networks need 4.6% fewer evaluations, with no change in outcome.

- A compressible `EjectorElement` no longer converges to a spurious zero-flow state with its mixing-plane static at ~1e13 Pa. The stress sweep found it in 9 of 11 "converged" cold draws, all reported as success. The C-D nozzle closure behind R0/R1 is now strictly monotone in `P_py`:
  - above the supply it passes reverse flow, where it used to return zero flow flat in `P_py`;
  - below the critical ratio the exit Mach is held at 1, where the supersonic-exit flux used to dip under the choked cap.

  Converged suite states moved by at most 4e-7 in primary flow (#619).

- A choked compressible `EjectorElement` now solves the subcritical regime (outlet between the critical back pressure P_c* and the dead-head P_b0). Before, it kept the critical entrainment there and pinned the outlet port at P_c*, below its boundary, so there was no root.
  - **Entrainment:** it now follows the Phase-1A Tier-1 droop, linear from omega_crit at P_c* to 0 at P_b0, and the outlet floats. The outlet pin switches in at P_b0.
  - **Diagnostics:** `diagnostics()` reports `p_b0_pa`.
  - **Critical-mode results:** unchanged.
  - **Cold stress draws:** 22 of 22 compressible draws converge, against 8 before.

- **A channel reports the friction factor its pressure drop uses.**
  - **Before:** every channel carries a convective surface, so its diagnostics reported the heat-transfer correlation's f at that correlation's own Re (16% apart at Re 70). That value kept the old Re floor sqrt(Re^2 + 64^2), so a stagnant duct showed f = 1 while its row used Poiseuille 64/Re with #634's floor.
  - **Now:** `f` is the residual's value: the C++ laminar-turbulent blend at the inlet Re, times the f multiplier.
  - **Enhanced surfaces** (ribbed, pin-fin), whose residual uses the surface f, still report it.
  - **The correlation's own value** is reported as `f_surface`.

- GUI: the hint a failed solve gives for a declared-direction tee (#640) now names what the Tee Inspector shows, such as the "Ideal momentum (parameter-free)" junction model and Flow Direction "Merge". It used to name solver classes (`IdealMomentumTeeElement`) and attributes (`flow_direction='merge'`), which nothing on the canvas maps back to.

- **A discharge-hole orifice between pipes loses Idelchik's in-pipe loss (#645).** This covers `IdelchikThick`/`Beveled`/`Rounded`, `McGreehanSchotsch` and `Lichtarowicz` between pipes. Before, their plenum-measured Cd took the velocity-of-approach factor E and the whole jet head was lost, so the flow came out 2-7% low at beta 0.3, 16-20% low at 0.5 and 34-38% low at 0.7.
  - **What changed:** the Cd is carried into the pipe as `Cd_wall * R(beta)`, and the jet now recovers into the far pipe like a metering plate's.
  - **Fidelity:** the loss matches Idelchik's diagrams 4-14..4-16 within 4.2% for beta <= 0.89.
  - **Unchanged:** plenum-fed holes, and the plate elements.
- GUI: the orifice Inspector shows a discharge hole's Cd factor in its pipe when it differs from 1.

- **A wall-heated compressible channel carries its heat in the momentum wherever it discharges (#672).**
  - #551 did this only into a plenum or chamber fed by the duct alone. Into a node that mixes other streams or heat, or into a boundary, the duct stayed adiabatic Fanno flow and the "Heating not in the momentum" flag was raised.
  - The solver now forms the duct's own exit state from its stream and its own wall heat, through the same mixer and Jacobian relay as a node, and the Fanno-Rayleigh kernel uses it.
  - The answer equals the explicit-outlet-plenum workaround to 1e-6. The flag remains only as a fallback, where no flow enters the receiving node.
- **A staged cold solve no longer switches the heat back on under its own adiabatic stage.** A solve nested inside stage 1 (the wall-free seed) cleared every duct's suspension on its way out. The staged solve now restores each duct's own flag, and nothing is staged with thermal coupling off.

- **Combustion composition derivatives are analytic (#716).**
  - **Complete combustion:** the products' Jacobian comes from the same code as their value (`complete_combustion_products_and_jacobian`), and T_ad's from the enthalpy balance. The forward difference inside C++ (eps 1e-6) read dT/dm_fuel 123% off at the phi = 1 softmin.
  - **Equilibrium:** implicit in each stage's extent and in the adiabatic root (`combustion_equilibrium_and_jacobian`). The forward differences (eps 1e-4 in Y, 2 K in T) were 4e-4 to 1e-3 off, and up to 12% in the product CO.
  - **Measured:** every combustor derivative within 7e-7 of central differences, phi 0.6-1.6.
- **The reforming equilibrium reaches its root near complete reforming (#716).** Its Newton took a finite-difference slope with an absolute step of 1e-8 mol, larger than the hydrocarbon left, and stalled off the root (a residual of -35 in ln Q at phi 0.9). It now uses the analytic slope inside a bracket. Only sub-ppm species move; T_ad is unchanged to 0.01 K.
- **Smoothed complete combustion conserves atoms inside its phi = 1 band (#716).** The O2 left is X_O2 - f n_req at every r, and the softmin applies on both sides of r = 1. Values move by at most about ln(2)/k (3.5e-5 of the fuel), within 1e-3 of phi = 1.

- The wall-coupling Jacobian now carries each wall side's sensitivity to its own composition. h, T_aw and the capacity rate |m| cp all move with the products behind a combustor, and the relay dropped that sensitivity, so a walled combustor's Jacobian was off by up to 18%. `ChannelResult` gains `dh_dY` and `dT_aw_dY` (filled by smooth channels). New `duct_htc_static_state_dY` and `cp_secant_mass_dY`, plus C++ `transport_and_dX` and `mole_gradient_to_mass` (#717).

- **A reversed area change is fed from the node the flow comes from (#731).** It read `from_node`'s state in both directions.
- **Area change composition derivatives are analytic (#731).** They were finite differences inside C++: gamma, R and the Wilke viscosity are now chained to the mass fractions.

- A merge chamber's Jacobian now carries both compositions. `chamber_merge_face_state` had no partials in the chamber gas's or the main stream's composition, so a side stream of another gas (fuel or O2 into air) left the main inlet's row up to 8.8e-3 off. `MergeFaceState` gains `dPf_dY`, `dPf_dY_main`, `dPtf_dY` and `dPtf_dY_main`. The momentum chamber closure's composition partials are analytic (they were a central difference), and a wall on a momentum chamber carries its static T's composition term (`MomentumChamberResult.dT_static_dY`) (#743).

- A GUI constant-K tee whose pair creates total pressure at its design branch share (`(1 - q) K_straight + q K_branch < 0`) is now refused when the network is built, with the cause and the fixes named, and the Inspector shows the same warning. Before, the solve converged and was then rejected by the second-law check with no hint why. Bassett's default joining pair does this at equal areas and theta <= 60 deg (both K = -0.04 at 30 deg, q = 0.5; measured about +0.05); typed K values are checked the same way. New `ConstantKTeeElement.design_net_loss` (#779).

- **Can-annular Bloch modes: wrong sign on the annulus admittance.** `can_annular_eigenmodes` summed the can admittance with the annulus admittance of the opposite sign convention, so its frequencies matched no physical system.
  - **Before:** for 24 cans (the docstring example) the push-push (m = 0) modes came out at 646.8 and 1198.6 Hz against 395.4, 831.5 and 1295.3 Hz from the exact symmetry reduction and a 1-D finite-element model of the full N-can graph (+64 %, one mode missed); push-pull was +8 %.
  - **Now:** every mode agrees with the finite-element reference to its discretisation error (N = 1, 5, 8, 12, 24; closed and open can tops). The dispersion relation is evaluated pole-free, which also recovers the modes with a pressure node at every junction (sin(k s) = 0) that the bare admittance sum cannot see, and the `|D|` acceptance test is now relative to the branch admittances as its header always said (an absolute 0.1 accepted any point).
  - `BlochMode.symmetry_type()` no longer labels m = (N-1)/2 "Push-Pull" for an odd can count.
- **`quarter_wave_resonator_tmm` ignored `A_tube`.** The branch admittance added a neck impedance normalised on the duct panel to a tube impedance normalised on the tube without converting either, so `A_tube` had no effect and `T21` was off by the porosity factor. It is now `1 / (rho c (z_neck / A_duct + z_tube / A_tube))`.
- **Orifice inertia used one face's end correction.** `orifice_impedance_with_flow` (and the liner helpers built on it) used `t + 0.425 d`; a perforate radiates from both faces, `t + 0.85 d`, the same 0.85 d `helmholtz_frequency` applies to a neck. The constant is now `kEndCorrectionFlangedBothFaces` in `acoustics.h`. It also refuses a non-positive or non-finite `c` or `rho` instead of returning garbage.
- **`tube_Q` was half the Kirchhoff value.** It returned `D / (4 delta_eff)`; the header and the wide-tube result give `D / (2 delta_eff)`.
- **`helmholtz_Q` had no physical basis.** It carried a spurious `sqrt(V / (A L_eff))` factor, so at a fixed frequency Q grew as the square root of the cavity volume (1.9x too high on the C++ test resonator). It is now the lumped `omega M / R = (d/2) L_eff / (delta_eff L_neck)`.

- **Annular-duct modes: the first azimuthal mode and the missed roots.** `annular_duct_modes_analytical` and `annular_duct_eigenmodes` (thin-annulus duct, rigid at x = 0) are now checked mode by mode against the closed form on 96 geometries.
  - **Missing n = 0 mode:** between rigid ends every m >= 1 has an axially uniform, purely azimuthal mode at f = m c / (pi D_mean), the classic first azimuthal mode of an annular combustor. Both functions left it out; it is now returned with `n_axial = 0`.
  - **Missed roots:** `annular_duct_eigenmodes` refined one Muller start per 50 Hz box and counted zeros minus poles as zeros, so it dropped any root sharing a box with a pole or another root (65 of the 96 geometries differed). It now adds the closed-form poles back, widens a box edge that sits on a pole, seeds Muller from the |D| minima, bisects contour steps that turn by more than pi/4, and evaluates the admittance by its series at cut-on (where it returned a 1e12 sentinel at the root). `n_axial` comes from the root, not from the order it was found in.
  - **`bc_ends` documented as it behaves:** the end x = 0 is always rigid, so `Open` gives closed-open (quarter-wave) modes, not open-open.
  - **f_max** is now exclusive in `annular_duct_modes_analytical` as in the numerical solver, the 200-mode cap is gone, and NaN or infinite `f_max`, `c`, `L` or diameters raise instead of returning an empty list.

- `AreaChangeElement` diagnostics now report the loss the network solved:
  `zeta` and `dP_static` no longer apply the kernels' `k1 = 1 + 0.35 M^2`,
  which the residual never uses (they read 3.6-5.0% high at a small-section
  Mach of 0.33-0.39). The reference `velocity` and `Re` are now taken at the
  smaller section, as labelled; a contraction used to report them at its
  larger inlet face (28 m/s against 114 m/s).

- An element discharging into a `PressureBoundary` reported the boundary's own set temperature and composition as its exit state (`Tt_out`, `T_out`), and a `mach_out` computed from that gas. A 600 K stream into an outlet set to 300 K reported 300 K. It now reports the gas that arrives, from the element's own stream and wall heat, including a reversed element discharging into its `from_node`. The solved pressures and flows were never affected (they do not depend on an outflow boundary's Tt), and the boundary's own result is still its setting.

- **Can-annular Bloch modes.** `can_annular_eigenmodes` is now checked root by root against an independently derived junction balance.
  - **Sign of the plenum admittance:** the annulus term was `(cos ks - cos psi)`; volume-flow balance gives `(cos psi - cos ks)`, so that for m = 0 the two half-sector stubs add compliance to the closed can instead of removing it. On the 12-can test case the old m = 0 modes were 60.2, 356.1, 735.0, 965.2 Hz against 315.6 and 763.0 Hz. The dispersion relation is now returned in pole-free product form, which also keeps the modes with a pressure node at every junction.
  - **Spurious modes:** a candidate was accepted when the dimensional |D| was below 0.1, which every point satisfies (A / (rho c) ~ 1e-5); one 12-can case returned 114 modes. Acceptance is now 1e-6 relative to the branch admittances.
  - **C++ `solve_argument_principle`** also seeds Muller from the |D| minima of each box, bisects contour steps that turn by more than pi/4, and accepts only points that are roots; it dropped 4 of 21 modes and returned non-roots before. `refine_root_muller` starts from a stencil of 0.1% of the guess instead of 1%, which jumped across a pole to the neighbouring root.

- **Canvas cards show inherited geometry (#559).** An area change whose faces inherit always read "F₁/F₀: 2.00" (placeholder areas 0.01/0.02) whatever the real ratio; it now shows the solver's resolved ratio. A channel that inherits its D showed a blank "D:"; it now shows the inherited D greyed (or an amber "?" when nothing can be inherited).

- **`CombustionState.fuel_burn_fraction` is the share of the fuel burnt.** Complete combustion burns 1/phi of the fuel when rich and leaves the rest in the products, but the field read 1.0 at every phi. It now reads min(1, 1/phi) (0.5 at phi = 2).
- **`combustion_state_from_streams` puts the products at the mixed pressure.** The reactant temperature already came from `mix()` at the lowest inlet pressure, but the reactants and products were placed at the fuel stream's pressure. Both now use the `mix()` pressure, so a 3 bar fuel line feeding 1 bar air gives 1 bar products.

- **Combustion helpers refuse a NaN instead of returning a wrong answer.** A NaN passes every `x <= bound` guard, so it slipped through.
  - **Before:** `set_fuel_stream_for_Tad(nan, ...)` returned the stoichiometric fuel flow, and the other nine `set_{fuel,oxidizer}_stream_for_*` solvers did the same. `fuel_lhv_molar(X, nan)` returned 0.0. `set_equivalence_ratio_mole/_mass(nan or inf)`, `complete_combustion_to_CO2_H2O`, `complete_combustion_products_and_jacobian`, `equivalence_ratio` and `oxygen_required_per_*_mixture` with a NaN fraction returned NaN. `complete_combustion(_isothermal)` accepted T <= 0 and P <= 0: a -10 K reactant gave T_ad = 1235 K.
  - **Now:** each raises `ValueError` (`std::invalid_argument`). `fuel_lhv_*` also refuses negative fractions, and `combustion_state_from_streams` refuses an oxidizer flow <= 0 or a negative fuel flow.

- A `CombustorNode` with an `EnergyBoundary(fraction=...)` now reports the heat that boundary takes as `{node}.Q_fraction`: the products' sensible heat at T_ad, times the fraction. It was never written for a combustor. The network ledger then saw an unexplained imbalance and demoted a correct root ("network balance does not close"), so any fraction loss on a combustor made the solve fail.
- The combustor's Jacobian with a post-combustion heat shift (wall heat or `fraction`) now includes the shift's dependence on the burnt composition at fixed enthalpy. A combustor fed by separate fuel and air streams had dT/dm off by up to 8% in a network with heat applied. The composition finite differences also now difference against the adiabatic temperature rather than the shifted one.

- **`combaero.network.diagnostics` imports without matplotlib.** matplotlib is the optional `examples` extra, but the module imported `pyplot` at load time, so a plain `pip install combaero` could not import it (found by the 0.7.0 clean-venv wheel check). Only the heatmap functions need it now, and they say how to install it.

- A compressible `EjectorElement` built through the Python API now gets the physics warm start the GUI gave it, so its cold solves converge like GUI-built ones. The warm start covers the port chamber pressures, `P_jct` (= P_c*) and the port-chain flows, taken from the closed forms.
  - **Where it lives:** `EjectorElement.seed_warm_start`, run at every `resolve_topology`. It moved there from `gui/backend/graph_builder`. It seeds only where no `initial_guess` was given, replaces its own seeds on a re-solve, and seeds nothing for the incompressible regime.
  - **Effect:** cold compressible stress draws went from 2 of 22 converged to 8 of 22, every critical one.

- **Ejector warm-start now reaches a boundary through a pipe on a port.** A channel on an ejector port makes the port an auto-inserted node with no canvas edges, so the warm-start found no boundary and seeded nothing there. The unchoked jet-pump case with an outlet pipe now converges (it stalled before). The double-choked case behind an outlet pipe still often fails to converge; that predates this change and is tracked in #567.

- **An ejector's auto-inserted outlet port chamber is sized at A_3 (#567).** When both the mixing area and the outlet pipe's D were left to inherit, the port chamber's area search hopped across the ejector to the primary inlet chamber (0.05 m^2 in the GUI default). The outlet then ran as a near-plenum. It now takes A_3. A port sized by an explicit outlet pipe keeps the pipe's area.

- **A double-choked ejector behind an outlet pipe now converges, and to the right root (#567).** The ejector warm-start seeded port pressures but left every flow on the port chains at the generic cold default -- tens of kg/s for a ~0.1 kg/s ejector -- which sent Newton to Mach-20 port states. It now seeds each port chain's flows from the closed-form critical operating point (mp, omega*mp, (1+omega)*mp), runs once the whole network is built (the outlet pipe may come after the ejector), and uses the outlet pipe's area when A_3 is inherited (the P_jct seed silently failed there). On the Huang T3-EH reference with an outlet pipe of 1.2-20x A_3, cold solves go from 6/16 to 12/16 converged, and the cases that previously "converged" to a wrong entrainment ratio (omega 0.29 / 1.73 instead of 0.79) now reach the critical root. The jet-pump regime is unchanged (15/16, no flow seed there). An outlet pipe close to A_3 itself still fails: that is a separate momentum-chamber closure limit, tracked in #567.

- A crossflow segment's reservoir entry (`entry_K`, `cb.channel_entry_drop`) applies to inflow only. Reversed, the duct discharges into the reservoir at static continuity, the segment's own convention at an exit. The mirrored `m|m|` form had charged the outflow `(1 + K)` dynamic heads of static drop, under-predicting any reversed flow through an effusion liner's duct entry (#471) or an impingement bypass inlet (#705). C1 through zero flow, so the Jacobian stays analytic.

- A starved or reversed impingement jet row is now crossflow-controlled:
  CR-3936's `eta_r` is taken at the crossflow over the jet's magnitude, so it
  tends to 1 from either side of zero jet flow (flagged `eta_r_extrapolated`).
  A reversed row used to drop to `eta_r = 0` (the reference temperature fell to
  the jet temperature), unflagged. Against CR-3630 Table D.2 end to end, the
  bypass `eta` RMSE falls from 0.163 to 0.078 (#467).

- **`fanno_max_length` returns the length to M = 1.** It bisected `fanno_channel` on "choked or exit M > 0.99", so it returned the length to M ~ 0.99: short of the closed-form L* by 0.8% at inlet M = 0.9 and 24.5% at M = 0.98 (argon, calorically perfect).
  - **Now:** the `fanno_length_between` Mach quadrature from the inlet Mach to 1 at the inlet's mass flux and stagnation temperature (the same constant-`f`, adiabatic, thermally perfect flow), matching the closed form to 1e-7 or better and the network's Fanno kernel exactly. `tol` and `max_iter` are unused.
- **`solve_P0_from_mdot` reaches small flows.** The bisection's lower end was `1.001 * P_back`, above the root of any flow smaller than a 0.1% drive passes, and such a target silently returned `1.001 * P_back` (+116% flow at 1e-3 kg/s through 2 cm^2 at 1 bar). The bracket now starts at `P_back`, and a target no P0 up to 1e14 `P_back` reaches raises instead of returning the last bracket.
- **`solve_P_back_from_mdot` converges on the flow, not only on P.** The bisection stopped once the bracket was `tol * P` wide, which at small flows (G ~ sqrt(P0 - P)) left the flow 73% off at 1e-4 of the choked flux and 1e-4 off at 1e-2 of it. It now stops on `tol * (P0 - P)`.
- **`bulk_velocity` and `stagnation_from_static` are exported from `combaero`.** Both were bound in `_core` and documented as `cb.bulk_velocity` / `cb.stagnation_from_static`, but `import combaero` did not expose them.
- **Incompressible `channel_flow` profiles carry `h` in J/kg.** `IncompressibleStation.h` is documented in J/kg (as `FannoStation.h`) but held the molar enthalpy, J/mol. The `channel_flow_rough` overload with a `k_loss_fn` now also puts the K-loss into its outlet pressure and profile, which were left at the friction drop alone (C++ API).
- **Units for the nozzle inverse solvers and `solve_orifice_mdot`.** `units_data.h` listed placeholders for `solve_A_eff_from_mdot` / `solve_P_back_from_mdot` / `solve_P0_from_mdot`, and `Pa` as the output of `solve_orifice_mdot`, which returns kg/s.

- **GUI backend: requests the builder could not represent are refused, not dropped or answered with a bare 500.**
  - **A node of unknown kind is refused.** It was skipped: a parallel branch of an unknown kind (for example from a file saved by another version) solved without that branch and reported success.
  - **An edge naming a node that is not on the canvas is refused,** naming the edge. It was ignored, or blamed on its element as "not fully connected".
  - **A flow edge drawn to a probe's handle no longer breaks the element** ("has multiple downstream links"); probes are display-only.
  - **A removed surface type ("dimpled") says it was removed (#339) and names the node,** as `map_surface_model` intended. The schema's surface union refused it first, with five pydantic literal errors.
  - **Invalid node or wall data names the node or wall and the field,** without pydantic's model name and documentation links.
  - **`/resolve_preview` writes a non-finite value as null,** as `/solve` does, instead of failing in the JSON encoder after the route returned.
  - **`/favicon.svg` and `/icons.svg` are a 404 without a built frontend,** not a 500.

- **GUI tells you when the solver backend is unreachable.**
  - Before: if the backend was not running, a solve attempt showed only "Solve failed. Check console for details." -- indistinguishable from a real solver failure.
  - Now: after a failed request the GUI probes `/health`; when the backend is down it says so and how to start it (`./scripts/run-gui.sh` or `combaero-gui`), instead of the generic message. A backend that is up but returns an error still shows that error's detail.

- GUI CSV export / `NetworkResult.to_dataframe()`: every row now carries the solve's status (`success`, `converged`, `consistent`, `outcome`, `final_norm`). Before, a timed-out or second-law-rejected solve exported numbers that looked valid, with `success` on element rows only. The GUI also warns when the export's own solve failed, and starting an export no longer clears the displayed results.
- A pressure boundary's exported (and probed) static `P` is its `Pt`: it is a reservoir. It used to be a made-up 101325 Pa shown beside the real `Pt`.
- An impingement array's thermal wall is exported once, as the aggregate that keeps the per-row values in `rows_Q` / `rows_T_hot`. The per-row walls were exported too, so summing `Q` over the export's wall rows gave twice the array's heat.
- A failed export shows the backend's message ("Export failed: ...") instead of "Check console for details".

- **GUI: node text stays inside its node.**
  - **Before:** the Effusion and Orifice nodes (120 px) and the Impingement node (150 px) printed their result rows past the node's end, which is most visible when the node is rotated.
  - **Now:** those nodes are 140, 140 and 172 px wide. On every node, a row that is still too long is cut with an ellipsis instead of overflowing.

- GUI: an orifice's bevel depth (Idelchik 4-18b) now reaches the solver. The builder dropped it, so a beveled hole drawn in the GUI was always solved at zero bevel depth.

- **GUI Probe and Live Telemetry: what they show is current, and every key has a unit.**
  - **A failed solve request no longer leaves the previous solution on the canvas looking current.** Every node and wall kept its old result after a request error (backend down, timeout), and the Probe, the Live Telemetry and the node cards showed it; a node a solve gave no result kept its old one too.
  - **A solve that did not converge is flagged** in the Probe ("not converged") and in the Live Telemetry: the backend returns its last iterate with `success: false`, which read as an answer.
  - **Picking a probe's target no longer clears every result in the network.** Drawing a connection to a probe now sets its target; it used to add a flow edge that changed nothing the probe showed and cleared all results.
  - **A Display Field you uncheck stays unchecked at the next solve;** every solve re-added it.
  - **Before the first solve the Probe offers the keys a solve gives the target's kind,** from the backend (`resultKeys.generated.json`, written by `scripts/gui_result_keys.py`), instead of a fixed list of node-state keys most elements never report.
  - **A result key with no catalogue entry shows its SI unit,** in the Probe, the Live Telemetry, the node cards and the Display Fields panel (about 130 keys, among them every `port_<i>_*` key of the tee and the ejector, the ejector's entrainment ratio and the thermal wall's UA, NTU and effectiveness, showed a bare number).
  - **The ejector's `port_<i>_*` keys carry their arm (P/S/O),** as the tee's carry C/S/B. Both port orders are now checked against a solve (`test_gui_port_arms.py`).
  - **The area-change card reads the solved `area_ratio`;** it read `ratio`, a key renamed earlier, and always fell back to the preview.
  - **The CSV export labels molar mass `mw [g/mol]`;** it said kg/mol (`State::mw` is g/mol).

- GUI tee inspector: the Constant K description now shows only for the Constant K model, and the straight-arm note says every junction model takes F_S = F_C instead of attributing it to Mynard alone.
- GUI numeric fields no longer push a value above their maximum to the node while typing (the tee's branch label briefly read "900°"); the value is still clamped on blur as before.

- **Humid air: `wet_bulb_temperature` and `dewpoint` no longer fail or lie at the edges of their range.**
  - **`wet_bulb_temperature` converges when T is above the boiling point at P.** Newton started at Twb = T, where the saturation humidity ratio is undefined once P_ws(T) >= P, and raised "did not converge" (20 of 1110 grid points, e.g. 375 K at 1 atm, 355 K at 0.5 bar). It now keeps a [lo, hi] bracket with bisection fallback; every previously converged value is unchanged and the new ones are within 0.55 K of CoolProp HAPropsSI.
  - **`dewpoint` refuses a dewpoint below the saturation-pressure range (IndexError)** instead of returning the -100 C clamp as converged (175 K at 10% RH returned 173.15 K, where P_ws is 6.9x the vapour pressure); between 163.15 and 173.15 K it now solves.
- **Documentation:** `humid_air_enthalpy_nasa9` is on the NASA-9 datum, not the psychrometric one (the two differ by hundreds of kJ/kg); `saturation_vapor_pressure` raises IndexError, not RuntimeError, and is over ice below 0 C; `docs/API_CPP.md`'s `HumidAir` listing matches the struct; `units_data.h` labels the psychrometric enthalpies J/kg dry air.

- **Inheritance labels from the closing audit (#559).** A discrete loss that takes its area from the upstream channel now says so (it showed no source); a channel's unset `Dh` is reported as the circular `Dh = D` rather than "user". Labels only; resolved values are unchanged.

- **The Inspector shows the geometry the solver will actually use (#559).**
  - Inherited diameters and areas now come from the solver's own resolve pass (`/resolve_preview`) instead of five separate frontend graph walks, so the placeholder, the "inherited from" note and the value used in the solve are the same number.
  - Fixes the mismatches those walks introduced: no fake `0.1` / `0.01` / `0.02` placeholders where the solver would refuse; an inherited channel `Dh` is now circular (= `D`), as the solver uses it, rather than an upstream channel's non-circular `Dh` it ignores; tee areas no longer "inherit" through plenums or multi-hop channel chains the builder does not follow; combustor and momentum-chamber areas show their inherited value instead of a fixed `0.1`.
  - Geometry the solver cannot resolve is flagged in amber with what to set, before Solve; defaulted values (chamber 0.1 m², tee `F_C` 0.01 m²) are marked as defaults; an inherited value is shown as a greyed placeholder in the current unit (like `Dh`'s "Circular (Auto)") rather than written into the field at full precision, where it read as user input; typing overrides it.
  - A newly dropped combustor or momentum chamber no longer stores `area: 0.1`, which counted as a user value and blocked inheritance (the solver used 0.1 too); the area now inherits like `Dh`. Existing saved networks keep their stored 0.1 until you click "Reset to inherited". The chambers also report where an area derived from `Dh` came from (`channel 'x'` / `derived from Dh`), which was previously blank.

- The isentropic, energy and stagnation Newton solves (orifice and nozzle
  kernels, `T0_from_static`, `T_from_stagnation`) no longer cap T at 6000 K:
  an orifice fed by a 6872 K plenum read its density at 6000 K and passed 3.6x
  its flow. The stagnation conversions' 200 K floor is `thermo::T_FLOOR` (10 K),
  as the others', and so are the combustion closures' 134 K inlet floor and the
  quasi-1D nozzle's 150 K one.

- **ISO 5167 metering Cd: small pipes, and the network's 'Stolz'.**
  - **Small pipes:** Reader-Harris/Gallagher now adds ISO 5167-2:2003's small-pipe term, +0.011 (0.75 - beta)(2.8 - D/25.4), below D = 71.12 mm. A beta-0.5 plate in a 50 mm pipe reads Cd 0.0023 higher.
  - **The network's 'Stolz':** the network `OrificeElement` computed 'Stolz' as a four-term Reader-Harris/Gallagher typed out in Python, with no Reynolds term from Stolz and no tap terms. It now evaluates the real ISO 5167:1980 Stolz equation in C++, through `metering_cd_and_dRe`, so no metering formula is left in Python.

- **A junction with several inlets and several outlets feeds each outlet its
  own flow.** The solver gave every sink of a multi-source element each
  source's own flow. That is right for a merge, which has one outlet, but at a
  `MultiPortChamberElement` or `IdealMomentumTeeElement` with two or more
  inlets AND two or more outlets, each outlet chamber took its velocity head
  at the junction's whole throughput.
  - **Converged results were wrong:** in the stress `ideal_nport` draws, port
    flows moved by up to 2x on fully forward roots.
  - **Solves stalled:** a port flow crossing zero jumped its chamber closure
    by about 10 kPa, which stalled 13 of 40 draws. All 40 now converge.
  - **The GUI's 3-port tees** (one inlet or one outlet) were not affected.

- **A plenum or boundary wired straight to a tee or ejector port no longer fails to build.** The junction needs a momentum chamber on each port and refused anything else ("port ... must be a MomentumChamberNode"), although the builder's own message said ports must connect to plena or boundaries. Such a port now gets an auto-inserted momentum chamber joined to the node by a lossless link, exactly as a channel on the port gets; a momentum chamber you placed yourself is still used as-is. Found in the #559 live testing (boundary on an ejector's suction port).

- **A channel's laminar friction is Hagen-Poiseuille down to zero flow.**
  The network friction kernels (`haaland`, `serghides`, `colebrook`) blend
  64/Re into the turbulent fit around Re 2300.
  - **The old floor:** the laminar branch shared the turbulent fits' Re
    floor of 64, so f was 5% low at Re 200, 29% low at Re 64, and saturated at
    1 below that (64x low at Re 1). A channel's drop then went as m^2 at low
    flow instead of linearly, with zero slope at zero flow.
  - **`colebrook` had no laminar branch at all:** it was 63% low at Re 200
    and 54% high at Re 2000.
  - **Now:** the laminar branch has its own floor of Re 1e-6, and Colebrook
    blends like the other two.
  - **Unchanged:** the turbulent fits, to 1e-12 at Re >= 1e4. With the
    default `haaland`, f moves by more than 1% only below Re ~ 300.
  - **One visible side effect:** f reported at exactly zero flow is now
    64/1e-6 rather than 1. The drop there is still zero.

- **A pipe-fed orifice no longer counts the approach velocity twice.** With a pipe upstream, the orifice applies the velocity-of-approach factor E = 1/sqrt(1 - beta^4).
  - **Before:** it also drove from the upstream stagnation pressure. Behind a node that carries the pipe's velocity head (a `MomentumChamberNode`), that counted the head twice: the flow was 1.2% high at beta 0.5 and 6.2% high at beta 0.7.
  - **Now:** with E applied, both kernels drive from the upstream static pressure, which is ISO 5167's tapping differential. Behind a `PlenumNode` (Pt = P) nothing changes, so the node type no longer decides the flow.
  - **API:** `orifice_compressible_residuals_and_jacobian` takes `P_static_up` after `P_total_up`. `orifice_residuals_and_jacobian` takes an optional `P_rho_ref`, the density reference, whose column is `d_mdot_dP_rho_ref`.

- An orifice whose discharge coefficient reads a Reynolds number (the ISO metering correlations and the discharge-hole family) now carries the coefficient's dependence on the upstream temperature, through the viscosity, in its solver Jacobian. Behind a heated plenum the row was 1e-3 off; solved results do not change.

- **A compressible orifice with a velocity of approach expands from the upstream static temperature (#645).** With E applied the nozzle starts from the upstream static state, but it paired the static pressure with the node's stagnation temperature. Behind a node with velocity (a `MomentumChamberNode`), the flow came out sqrt(T_s/Tt) low: 0.06% at pipe Mach 0.08, 0.3% at Mach 0.18. The static temperature now lies on the isentrope through the node's (Tt, Pt). Behind a `PlenumNode` (P = Pt) the flow is unchanged.

- **The orifice's McGreehan-Schotsch option no longer claims a crossflow.** A plain `OrificeElement` is plenum-fed, so the correlation's inlet-crossflow term always sees U1/Vi = 0. The GUI option read "Cooling Hole with Crossflow" and now reads "plenum-fed", like the plate elements. Only `EffusionPlateElement`'s duct-fed panels feed it a crossflow.

- **The metering orifice's Jacobian carries its Cd's flow dependence.**
  `ReaderHarrisGallagher` (the `OrificeElement` default), `Stolz` and
  `Miller` read the pipe Reynolds number Re_D = 4|m|/(pi D mu), so the row
  `m - Cd(Re_D(m)) K` has slope `1 - K dCd/dm`. The Jacobian gave 1: about
  1% off for a 2 cm orifice at 1600 K, in both regimes. The documented gap is
  closed: the C++ `Cd_sharp_thin_plate_and_dRe` (new, also in Python) gives
  the Reader-Harris/Gallagher derivative analytically, and Stolz and Miller
  get their closed-form derivatives. Results do not change; convergence
  paths can.

- Cooling correlation validity ranges given as a single tested value (Han, Zhang and Lee 1991 e/D 0.0625, P/e 10, W/H 1 and angle; Rallabandi 2009 W/H 1 and 45 deg; Taslim and Spring 1987 geometry) now flag values above that value as extrapolated. Previously `RibRange`/`ImpingementRange` read `hi == lo` as an open upper bound, so a set never flagged a value above its tested one. Bounds are now read with a 0.5% relative tolerance (`RibRange::rel_tol`) so geometry quoted to three figures (W/H 0.286 against 1/3.5) does not flag its own source's data.

- **`reforming_equilibrium` and `reforming_equilibrium_adiabatic` now reach the simultaneous equilibrium.** They reformed each hydrocarbon and then applied WGS once, so reforming was left short of equilibrium after WGS had moved CO, H2 and H2O.
  - **At 650-850 K, the pre-reformer range:** H2 was 7-49% low against Cantera. A CH4/steam feed at 900 K came out at 734.1 K instead of 711.9 K.
  - **Reforming extents could only run forward,** so methane never re-formed from CO + H2 (methanation). Ethane and propane at 600-800 K were left 0.2 off in mole fraction.
  - **The fix:** the reactions are re-equilibrated in rounds until nothing moves, with reforming extents bracketed in both directions. Over 600-1400 K and four feeds (CH4/steam, C2-C3 mixtures, steam-to-carbon 1.5 with CO2, syngas with CH4), composition is within 2.4e-4 of Cantera and adiabatic temperature within 0.12 K.
  - **`combustion_equilibrium`** uses the same solver on hot burned gas, where its results are unchanged for CH4. Its analytic derivatives (#716) iterate the same rounds.
  - **`converged` is now reported** instead of always True. Rich propane at phi >= 2.5 is steam-starved after complete combustion, the rounds do not converge, and the result is flagged rather than passed off as an equilibrium. That case was equally wrong before.

- GUI: the "Restored your last session" note above the canvas now goes away when a file is loaded or the canvas is cleared, instead of staying until it is dismissed by hand.

- **`NetworkRunner` solves a saved drawing the way the GUI does.** The GUI's solve and export routes turn on the Newton-homotopy last retry (#642), but `NetworkRunner.solve` and `.sweep` did not, so a script could fail on a network the canvas converged. All three now build their solver through `gui/backend/solver_setup.py`.

- A solve whose starting residual is non-finite no longer reports convergence. |F| = NaN failed every comparison in the best-iterate tracking and the success guard, so the attempt returned `__success__ = True` at a NaN point; now the start counts as unevaluable (outcome `ERROR`, `FloatingPointError` in strict mode, the rows' owners in `__caught_exceptions__`) and the retries take over. Seen on a compressible ejector retried from its incompressible proxy.

- **A network with no pressure difference now solves to zero flow.**
  - Before: two boundaries at equal total pressure (any network whose only solution is zero flow) failed to converge -- the default non-zero flow seed stalled on the flat `dP/dm` at `m=0` -- and the solver returned a wandering best iterate. In the GUI this surfaced as a bare "Solver: Failed".
  - Now: when the boundaries have no pressure driving potential, a failed cold solve is retried once from a zero-flow seed and accepted only if it genuinely converges, so the valid `m=0` solution is returned. A real pressure difference or an imposed mass flow is unaffected, and the retry never forces zero flow on a network that actually carries flow.

- The solver's stall detector catches a dead plateau just under its "near the root" gate (|F| < 100x tolerance). A cold attempt whose best |F| does not improve at all over three windows (15 s) now hands over, instead of grinding without end.
  - **Where it showed:** a 5-row hot impingement array's LM fallback sat at |F| = 0.0718 for 175 s with no timeout set.
  - **Crawls are untouched:** any improvement still counts as a crawl, and seeded attempts never run the detector.

- **An orifice whose Cd vanishes no longer raises.** A crossflow-fed hole's
  McGreehan-Schotsch Cd underflows to exactly 0 where its pressure drop
  vanishes, and the nozzle kernels refused `A_eff <= 0`, so the solver's
  safety net rejected those iterates (35 caught probes across the stress
  effusion family). Both kernels are exactly linear in Cd, so a Cd below
  1e-12 is now evaluated at Cd = 1 and scaled: a trace flow, bit for bit,
  instead of an exception.
- **A trace species no longer divides by zero in a pressure loss's
  Jacobian.** The composition column took a central step of half the mass
  fraction, which is 0 at a denormal fraction (5e-324). Below 2e-6 it now
  takes a forward step of 1e-6.

- A solve that converged with a compressible momentum chamber on the supersonic root of its closure (static pressure below the sonic pressure of its flow) is now reported INCONSISTENT and named in `__supersonic_chambers__`. It used to be reported as converged and consistent. Seen in the GUI as a channel into a narrow orifice: the junction between them at Mach 1.46, the channel "leaving" at Mach 1.67 unchoked, and the orifice's static pressure rising across it. The auto-retry then finds the subsonic root, with the same flow and total pressure. Two Wang 2014 junction endpoint points (area ratio 1.56, q = 1) had been scored on such a root and are now unscored (cross-source MAE 0.127 to 0.123 over 188 points; the floors were not tightened).

- GUI: the tee branch angle is now refused by the backend outside 0 < |theta| <= 90 deg, the rule the frontend already applied. A saved or API-supplied graph with 900 deg solved silently as a 180 deg tee, because the closures read the angle only through cos/sin; 0 deg (the branch lying on the straight run) is now refused on both sides, and the Inspector's help text no longer claims 0-180 deg.

- GUI: a node's `T` and an element's `T_in`/`T_out` are labelled "(total)". They hold the total temperature at every node, but `T_in`/`T_out` were listed under a "Static" heading. A momentum chamber's `T_static` is now labelled "T (static)". The key convention itself (T static, Tt total, as for pressure) is #782.
- GUI: plenum, combustor and momentum-chamber nodes no longer offer a `T` initial guess. Temperature is never a solver unknown, so the guess did nothing.

- A `ThermalWall` with a side that has no convective surface (a `PlenumNode`, an orifice, or a loss element, combustor or momentum chamber with no surface area) now fails validation, naming the side. Before, it exchanged no heat and raised nothing: the wall was silently inert.
- GUI: the orifice and area-change nodes lose their thermal handles. A saved network that walls one of them (or a plenum, effusion plate or film holes) loads with that wall removed and a note naming it. An effusion plate stays off the list on purpose: it owns its wall, and its bore is already in Andrews' internal correlation (the throat term on pi d L_hole).

- A thermal wall on a node (`CombustorNode`, `MomentumChamberNode`) now has an exact solver Jacobian. The wall relay put the node's own-flow heat-transfer derivative on the node's static-pressure column instead of its inflows, and the chamber's static temperature was held fixed in it.
- `CombustorNode`'s wall heat transfer is evaluated at the static state its throughput implies, as ducts are since #636. Reading the stagnation temperature as static put the adiabatic wall temperature above the stagnation temperature: +4.7 K at Mach 0.15, +46 K at Mach 0.49.
- A combustor carrying heat (an `EnergyBoundary` Q or a wall) now has the derivative of its outlet temperature with respect to its inflow: the heat is spread over the throughput, and that term was missing in both the complete and the equilibrium combustion paths.

- **The wall relay's pressure column (#496).** A wall's heat moves with its channel's static pressure at fixed flow through the adiabatic wall temperature: T_aw = T + r v^2/(2 cp) with v = mdot/(rho A), so dT_aw/dP = -2 (T_aw - T)/P. h does not move (Re, Pr and k do not depend on P).
  - **Before:** the relay had no pressure column. That was 1e-3 of the scaled Jacobian in a two-wall network, its largest entry after #463.
  - **Now:** `ChannelResult.dT_aw_dP` carries it through the smooth, ribbed, impingement and pin-fin paths, and the relay chains it into the feed node's P. The network Jacobian is exact to ~3e-9.

- **A merge chamber keeps its main-inlet face when the main flow reverses (#493).**
  - **Before:** the chamber hands its declared main inlet the face state of its impulse balance, but only while the main was an inflow. At m_main -> 0+ the face sits below the chamber pressure, because the side streams carry momentum out. So the main's exit pressure stepped 245 Pa, and its residual 3.5 kPa, as its flow crossed zero.
  - **Now:** the impulse holds for either sign (the x-momentum flux through the face is m_main^2/(rho A) both ways), so the face is kept. A reversed main gets the chamber's fluid at the face, and the outlet carries the side streams less the main's outflow.
  - **Measured:** continuous through zero; Jacobian exact (<= 5e-9 scaled) with the main forward, reversed and just past zero. Forward values are unchanged.

- **Mass-flow derivatives of a non-circular channel (#463).** `channel_smooth` took its mass flow on pi D^2/4, so for a channel with a separate hydraulic diameter dh/dmdot, ddP/dmdot and dT_aw/dmdot were off by (Dh/D)^2. Values were unaffected.
  - **Fix:** a new `flow_area` argument (default: the circular area) carries the element's real cross-section. `ConvectiveSurface` passes it on every path: smooth, ribbed, impingement and pin-fin.
  - **Measured:** global Jacobian error 1.2e-1 -> 1.2e-3 on a wall-coupled flat-duct network. The remainder is the missing pressure column, #496.
- **The wall relay through a feedback loop (#481).**
  - **Cause:** a wall coupling a node's outflow to its inflow feeds the node's temperature back into its own heat input, T' = a + b T', and the relay read that node's half-built sensitivities. Separately, a node read ahead of its order (a back-edge) gave an empty relay on every pass.
  - **Effect:** dT/d(cooling flow) was 19% short in a two-wall series network.
  - **Fix:** both reads now take the previous pass's sensitivities. The passes iterate until every borrowed sensitivity equals what the pass computed (1e-9).
  - **Cost:** junction networks read ahead through chains that carry nothing and settle in the first pass. The 240-case junction benchmark takes +5%; impingement and mixing are unchanged.

- **Jacobian completeness (#481 C1, C4, C5).** The analytic network Jacobian now equals differences where it dropped a dependence.
  - **Theta-sourced `PressureLossElement`:**
    - **Unburned temperature:** the combustor's unburned temperature enters through theta = T_b/T_u - 1 and as the correlation's reference temperature. It was never relayed and now is, through the combustor's inflows.
    - **Inlet composition:** a head loss's density depends on the inlet composition, which now has a composition column.
    - **Before:** errors were 2.6% in d/d(air flow) and 5% in d/d(fuel flow); both now hold to ~4e-6.
    - **Lazy:** the element's differenced Jacobian is now a `LazyJacobian`.
  - **Temperature clamp:** a node whose derived temperature sits on the [200, 5000] K clamp relays zero temperature sensitivity. Its flat value used to carry the mixer's full slope, 5% off.
    - **New check:** a converged solve with a flowing node on the clamp is now reported INCONSISTENT (`__success__` False, nodes in `__T_clamped__`). The clamp discards energy -- 190 kW in a heated test plenum -- so the residuals can close while the balance does not.
    - **Not flagged:** a stagnant node on the clamp, whose heat is already reported as `Q_withheld`.
  - **Wall heat relay:** this took dT/dQ = 1/(cp * sum m). The mixer spreads Q over m_eff: |m| above `MIXER_HEAT_MDOT_FLOOR`, a C1 parabola below. The old form was 2x off below the floor and had the wrong sign for a negative total. `MixerResult.dT_mix_dQ` now gives the exact value, also for the combustor mixers.

- **No throw, no NaN and no jump where a solver probe can reach (#481).** The suite's converged solves (485) were probed through flow reversal, zero flow and out-of-range pressures. Each fix is exact wherever its closure is defined.
  - **Turbulent-only Nusselt correlations** (`dittus_boelter`, `sieder_tate`, `petukhov`) threw for 2300 <= Re < 1e4 in `channel_smooth` and `htc_circular_channel`. They now use the VDI Heat Atlas transition: linear in Re from laminar Nu at 2300 to the correlation at 1e4, with exact derivatives. Values from 1e4 up are unchanged.
  - **`fanno_channel_flow`** threw at a back pressure <= 0. That now gives the choked flux, continuous with P -> 0+.
  - **Ejector:** the residual was NaN for every primary flow below ~0.9 of the operating one, and on most pressure probes. The jet-pump and critical closures took sqrt of a negative once the mixing-plane pressure passed a supply's, and a hard M_py = 0 switch divided by zero.
    - **Fix:** their arguments continue through a C1 positive floor, `dpos_floor`: eps^2/(2 eps - a) below eps. Unlike an exponential tail at the small eps a domain guard needs, it never underflows.
    - **Also:** the subsonic mixing root is now taken in its cancellation-free form.
  - **Legacy `TeeJunctionElement`** went NaN with an idle branch: the datum is the collector there (x = 1, phi = 0), and K -> 0 but evaluated 0 * inf.
  - **`MultiPortChamberElement`** classified a wrong-direction port by a 1e-9 kg/s tolerance but its suppliers and collectors by 1e-9 m/s. An inlet carrying 1e-9 kg/s outward passed as "not wrong", left no supplier, and the kernel raised. Both tests now use the velocity.
  - **A node with no inflow** (all its flows leaving) set its flow to 1 kg/s. That put a 4.3 kPa jump in a momentum chamber's dynamic-head row as its last feed reversed. It now carries no flow, continuous with a vanishing feed.
    - **Not shipped:** a continuous T/Y fallback for such a node, mixing the states it would be fed from weighted 1/|m|. It removed most remaining zero-flow jumps, but doubled the cold-start time of hot impingement arrays and failed one, with no gain elsewhere. Its T/Y stay at the 300 K-air fallback.
  - **Not changed:** a node's state still switches source when a series line reverses (upwinding). The merge chamber's main-inlet face still switches off for a reversed main; that is tracked as #493. In the 240-case junction sweep, none of the 10 solvable failures ends near zero flow.
  - **Measured:** over the suite's 490 converged solves, every throw and NaN is gone and zero-flow jump rows fall from ~750 to 543. On the 292-network benchmark the outcomes are unchanged, and junctions take +2% evaluations: a just-reversed junction port now takes the soft barrier, not the kernel's refusal.

- **Compressible channel derivatives are exact at low drive; the flow is smooth below the 1 Pa floor (#481).**
  - **Cause:** at low drive the duct's inlet and exit Mach nearly coincide, with (Me - M_in)/Me ~ drive/P (7e-6 at 1.5 Pa on 2 bar). `fanno_channel_flow`'s implicit derivatives combine terms ~P/drive times larger than their result, and two of those terms were finite differences: the inlet flux slope and the exit-pressure-per-flux slope. Their 1e-8 to 1e-10 rounding, amplified ~1e5, put up to 8e-3 into dG/dPt0 and dG/dP_target at a few Pa. Below the floor it also put 1e-5 noise into the flow itself, through the blend exponent sigma, so that G jumped when Tt or the composition was nudged.
  - **Fix:** both slopes are now closed form (energy plus isentropic relations); the moving inlet limit is taken by Leibniz's rule.
  - **Result:**
    - dG/dPt0, dG/dP_target and dG/dTt0 hold to <= 1e-5 from 1.5 Pa of drive up (were up to 8e-3, 8e-3 and 5e-3);
    - G's noise below the floor fell from 1e-5 to 1e-8 of G;
    - the kernel is 13% faster (4.75 -> 4.1 ms per call), and compressible-duct network solves ~10% faster.
  - **Tests:** reverting either slope to its difference fails the new tests.

- **Collision integrals are exact and smooth; transport derivatives are analytic (#485).**
  - **Interpolation error:** Omega*(2,2) was interpolated linearly in raw T* on the Monchick-Mason table. That was off by up to 2.7% between nodes and left 2-17% kinks in dmu/dT at every node.
  - **Table error:** the table's own Lennard-Jones column was off by up to 0.6% at the nodes; Kim & Monroe (2014) and an independent quadrature agree with each other to 5e-5.
  - **Now:** the delta* = 0 column (values and slopes) comes from quadrature of the classical collision integrals, generated into `include/collision_integral_data.h` (`thermo_data_generator/collision_integrals.py`). It agrees with Kim & Monroe to <= 1e-4 from T* = 0.3 to 400, between nodes included. The grid extends to T* = 400, which covers hydrogen to about 15000 K.
  - **Interpolant:** cubic Hermite in (ln T*, ln Omega*) with the tabulated slopes; C1 in T. Picked on a study against the exact function: 0.032% max between nodes, against 0.093% for the best spline and 0.23% for PCHIP.
  - **Polar columns:** scaled by the same exact/tabulated ratio, so they no longer sit 0.6% off the LJ column at T* = 100.
  - **Derivative:** `viscosity_and_jacobians` now returns the analytic dmu/dT instead of a +/-1 mK difference.
  - **Measured effect:** the old kinks cost no convergence (24,039 solves: evaluations -0.05%, no outcome flips). They did cost a 0.3% Jacobian error near zero drive in the compressible channel's dG/dTt0 (#486 had attributed it to the friction laminar blend), and a ~1e-8 precision floor on Re-dependent integrals.
  - **Values** change by tenths of a percent: air +0.1% at 300 K; 1 of 240 cooling validation series moves by 0.1 point.

- **Air's viscosity was 11.8% low, O2's 44.5% (#485).** The generated transport data took Lennard-Jones parameters from NUIGMech1.1 entries that its own authors flag as unreferenced (`theoret trans`). For O2 that meant eps 676 K, against 107.4 K in GRI30, JetSurf2, Aramco2 and San Diego. H2 (304.7 vs 38.0 K), H2O and Ar were wrong too.
  - **Generator change:** the extractor now lets a flagged transport set yield to the first later unflagged source, for transport only. The header is regenerated; thermo is unchanged.
  - **Viscosity at 300 K, after:** O2 -0.2%, N2 +1.0%, Ar +1.6%, H2 +0.4%, CO2 +0.3%, air +0.6% against literature. Alkanes move 1-3%.
  - **Affects every Re, Pr and Nu computed for air or O2-bearing gas.**
  - **Validation:** 7 of 240 cooling series move, all Andrews effusion, because they compute Re from rig conditions:
    - internal h accuracy bias -13.5% -> -9.9%;
    - fidelity pooled bias -10.4% -> -7.1%;
    - plate B -19.6% -> -16.5%;
    - overall eta plate B +23.7% -> +24.6%.

    The documents quoting them are updated. One note had blamed hand-typed properties for a gap the library's own viscosity caused; it is corrected.
  - **Tests:** the Cantera transport tolerance drops from 35% (which let this through) to 7%. A literature test pins pure-species and air viscosity independently of any mechanism. Restoring the old data fails 4 of 6 literature tests and 9 of 17 Cantera transport tests.

- **The compressible `ChannelElement` passes its choked flow, not more (#481).** Its residual is now in mass-flow form, `m - m_calc`. `m_calc` comes from the new `fanno_channel_flow`; it rises monotonically as the back pressure falls and saturates exactly at the choked flow, like the compressible orifice. Measured on a 20 mm by 1 m duct at 2 bar:
  - **Choked flow:** pressure-driven ducts past choke converged 5-7% (exit lost) and 6-20% (recovered) ABOVE their choked flow, because the infeasible region was patched with a barrier. They now converge on it exactly (1e-8).
  - **Failed solves:** the 60 kPa back-pressure case now converges.
  - **Unchoked flows:** unchanged, being the same Fanno physics.
  - **Imposed flow above choke:** it raises the supply pressure until the duct can carry it.
- **Compressible channel defects the redesign removes (#481):**
  - **Reversed flow** marched from the upstream node; it is now fed from the downstream node's stagnation state.
  - **The composition Jacobian** `d_dP_dY` was off by about 1e5x, because it differenced a different function. The new form has no composition column (see Not changed in #481).
  - **`get_spatial_profile`** raised TypeError for every compressible channel. It now marches from the duct's inlet static state.
  - **Residual shape past sonic:** dP(m) was non-monotone, jumped, and had d/dm = 0 past inlet-sonic.
- **Residual-row scaling is declared by the element (`residual_scale_kind`)** and documented once, under the searchable tag `RESIDUAL-ROW-SCALING` in `NetworkSolver._build_residual_scales`.
  - **The bug:** the compressible channel's new kg/s row, filed as a pressure row, stalled a GUI tee network (697 evaluations, failed). Correctly filed, it converges in 48 evaluations and 1.8 s; the old element took 1.9 s.
  - **Also fixed:** the scaling comment said `F_real * D_f`, the inverse of what the code does.

- **A `MassFlowBoundary` between elements now injects its own stream (#481).** The mass balance counted its `m_dot`, but the mixing dropped its enthalpy and species. Measured: 50 g/s injected at 800 K left the node at the inflow's 300 K (-26 kW).
- **`EnergyBoundary.fraction` now scales sensible enthalpy (#481).**
  - **What it means:** h(T) - h(298.15 K) at the stream's own composition; on a combustor, the products'.
  - **Before:** it scaled absolute enthalpy, which includes formation enthalpy, so a "-5%" HEATED air at 298 K and combustion products.
  - **Now:** -5% removes exactly 5% of the sensible heat. The heat applied is reported as `{node}.Q_fraction`.
- **Mixer Jacobian with heat or `fraction`.** Measured against finite differences, 1e-9 now:
  - d/dY used the base enthalpy where the final one belongs: 46% off with Q.
  - The fraction term multiplied the composition term by (1 + f): d/dm off by 2%, d/dY by 280%.
- **Heat is spread over the real flow down to 1 mg/s (#481).** The mixer floor was 2 g/s, so a node at 1 g/s took up only 80% of its heat, unreported.
  - `MIXER_HEAT_MDOT_FLOOR` is now 1e-6 kg/s; convergence measured the same.
  - What a stagnant node cannot take is reported as `{node}.Q_withheld`, as is wall heat aimed at a `WallNode` (previously booked as `Q_wall_out` with no stream to carry it).

- **Reversed flow conserves energy, and its temperatures are right (#481).** A node used to mix only the elements declared INTO it, at whatever sign their flow had.
  - **Upstream element reversed:** it entered as negative mass at the wrong node's temperature. Measured: feeds at 300 K and 800 K gave a 1380 K node (+57 kW).
  - **Downstream element reversed:** the fluid it delivered was never mixed in. Measured: a 600 K feed left the plenum at 400 K.
  - **Also wrong:** densities, and with them the flows.
  - **Now:** streams follow the sign of a 2-port element's flow, the propagation order is built from the signed flows, and walls put their heat into the node the flow goes into.
  - **Unchanged:** forward-flow networks are bit-identical.
- **The residual is a function of x alone (#481).** Two kinds of state used to carry over from the PREVIOUS evaluation, whatever trial point the solver had probed:
  - **Walls on back-edges** read ahead of the propagation order, so the same x evaluated three times differed by 3 Pa, then 0.003 Pa. That also corrupted finite-difference Jacobians (38% on a GUI tee network).
  - **k(T) wall layers.**
  - **Now:** the propagation repeats until such states settle, and k(T) is iterated at the call's own wall temperatures.
  - A mixing-plenum test marked "non-deterministic in suite" now passes reliably; its xfail is removed.
- **Incompressible channel Jacobian under reversed flow had the wrong sign** (+9.7e4 vs -8.0e4 by finite differences). dP is odd in m_dot, so its slope is even. The wall relay had the same error: surface correlations take |m_dot|, and their dh/dm was relayed unsigned (+0.096 vs -0.096).

- **A non-positive htc no longer stops the wall (#481).** A correlation past its range can return h <= 0 at an iterate; Florschuetz's jet-array bracket does for a nearly stopped jet under crossflow.
  - **Before:** `wall_temperature_profile` threw, and `wall_coupling_and_jacobian` switched to R = 1e15 while taking dQ/dh from the raw h.
  - **Now:** h is floored smoothly below `WALL_HTC_KNEE` = 1e-2 W/(m^2 K) by knee^2/(2 knee - h), which is C1, positive and never underflows. dQ/dh is the slope of the floored value.
- **Iterates the model refuses are rejected (#481).**
  - **Before:** an exception was answered with F = x - x_best and J = I. That F is SMALL near x_best, so the trust region could accept the step. NaN and inf went straight to MINPACK, and any recovered exception turned a failed solve into an ERROR.
  - **Now:** the answer is a residual ten times the last physical one, with its Jacobian, and NaN and inf count as refusals. A failed solve lists its rejected probes. A model that fails on every evaluation is still an ERROR.
- **With `thermal_coupling_enabled = False` the walls report no heat.** They reported Q and `Q_wall_out` for heat the solution never saw.
- **Cold convergence, impingement array under a hot duct** (4 row counts x 4 pressures, no retry): 9/16 after the htc floor and channel sign, 11/16 with the penalty, 16/16 with sign-aware streams. The failures sat in reversed-jet basins at 1500-3100 K, which the old mixing made self-consistent.
- **Every converged solve in the test suite is now checked for global energy and node mass closure** (`python/tests/conftest.py`). A test that builds a deliberately non-closing network opts out with `@pytest.mark.no_closure_check`.

- **Network energy balance now closes for every cooling configuration (#471).** Three leaks, each locked by `test_energy_conservation_cooling.py`:
  - **The mixers withheld part of the heat they were given.** They divided Q by `sqrt(m^2 + 1e-6)` instead of the flow, which lost 5e-5 of Q at 0.1 kg/s, 0.5% at 10 g/s and 11% at 2 g/s. This affected `mixer_from_streams_and_jacobians` and the two combustion variants. Q / m is now exact at or above 2 g/s (`kMixerHeatMdotFloor`); below that, a C1 floor applies.
  - **Wall heat put on a boundary disappeared.** A boundary's state is fixed, so that heat leaves with the stream crossing it. It is now reported as `{wall}.Q_to_boundary` and `{boundary}.Q_wall_out` (W). Before, the hot duct of an impingement array lost its walls' entire Q from the balance.
  - **Re-solving a network counted wall heat twice.** Every solver setup hung another wall `EnergyBoundary` on the heated nodes beside the previous one, which kept the last attempt's heat. That happened with a second `NetworkSolver` on the same network and with the solver's own retries. Setup now replaces it.
  - **Result:** closure is 1e-9 relative or better, against 1e-5 to 1 before.
- **A failed cold solve of a network with walls is retried from its wall-free flows.** The walls move a lot of heat into small flows. A cold Newton step that reversed one of them left a node with heat and nothing to carry it.
  - **Impingement array under a hot duct:** 6 of 9 cases converged cold, now 9 of 9 (3/5/8 rows, 1-5 kPa).
  - **Old mixer, for comparison:** 7 of 9. The exact heat made the cold path no easier.

- **The state relay of a single-source, multi-sink element took the stream's flow derivative at the wrong node.** This affects branching junctions. The stream into a sink is `flow_at_node(sink)`, but the relay differentiated `flow_at_node(source)`.
  - **Effect:** a heated node behind a branching tee's straight leg saw d/d(m_com) only and lost the −1 on m_branch. That Jacobian entry was 0 against a finite difference of −1.9.
  - **Now:** the global Jacobian matches finite differences to 2e-5. Solutions are unchanged; only Newton's path is.

- **The centred side-stream term dropped the segment's own-flow density term (#471).** The term is `m^2 (1/rho_from - 1/rho_to)/(2A^2)`. It is restored by building every segment from two half-stations at their nodes' densities.
  - **Florschuetz chain vs Eq. 8:** max 4.1% -> 3.8%.
  - **Replaced API:** `side_stream_momentum_drop` (unreleased) is replaced by `station_half_drop`.

- **`critical_pressure_ratio` returned 0.5636 whatever the gas or temperature (#471).** Its golden-section search kept its probe points in the opposite order to its update and shrank onto a fixed point.
  - **True value for air:** 0.528 at 300 K, 0.542 at 1500 K. That is where the isentropic mass flux peaks at M = 1.
  - **Choked mass flux:** it came out 0.1-0.3% low, and `nozzle_flow` declared choking early. The compressible orifice uses `nozzle_flow`.
  - **Tests:** `compressible.h` had no C++ tests; the new ones check M = 1 and the flux peak.
- **The compressible `ChannelElement` filed its inlet-pressure sensitivity under the upstream static P (#471).** The Fanno march inverts the inlet static state from the inlet TOTAL pressure (#359), so the C++ derivative is with respect to Pt.
  - **Behind a plenum:** harmless, since P and Pt are tied.
  - **Behind a momentum chamber:** 25% wrong, because the two are separate unknowns.
  - **Now:** the global Jacobian matches finite differences to 5e-6 there.

- **Orifice Jacobian now carries the discharge-hole Cd's flow dependence (#471).** For Idelchik, Lichtarowicz and McGreehan-Schotsch, Cd depends on the hole Reynolds number, hence on the flow. The residual's d/d(m_dot) left that out and was off by 1.5-3.6%. It is now analytic via `discharge_cd_and_derivatives`. Solutions are unchanged; Newton converges where it previously stalled, for example 30-degree effusion into a merge chamber.

- **Impingement channels in the GUI are a working example again (#462).**
  Since #460 the convective area counts the impingement holes, but channels
  defaulted it to the wetted wall `pi D L`. On a default 1 m channel that
  spread the flow over about 820 holes: `Re_j` 318, flagged.
  - **Default convective area by model** (`ChannelElement.default_convective_area`):
    - impingement row: footprint matched to the channel cross-section,
      `A_flow (xn/d)/(z/d)`;
    - single jet: Goldstein's averaging disc, `pi (R_D d)^2`;
    - everything else: unchanged.
  - **Result:** the GUI default array at 0.01 kg/s now gives `Re_j` 6,373,
    inside Florschuetz's 2,500-70,000.
  - **New channel diagnostics** `Re_surface` and `surface_extrapolated`,
    shown in the GUI as "Re (surface)" and "Surface out of range".
  - **Docstring correction:** the single-jet Nu was described as LOCAL at
    R/D. Han's Eq. 4.1 is an AVERAGE over the disc of radius R.

- **Impingement elements used the wrong jet flow (#460).** Both
  `ImpingementModel` and `SingleJetImpingementModel` took the jet flow as
  `rho v * ConvectiveSurface.area`: the convective area scaled the
  element's flow by `A_surface/A_cross`. Both sources define Re on the jet's
  own flow through its own hole (Goldstein 1986: nozzle Re; Florschuetz
  1981: jet mass velocity on the hole area).
  - **What changed:** the element's mass flow is now the jet flow. The
    convective area only counts the array's holes; for the single jet it is
    only the patch the reported h applies to.
  - **Validated against data.** Driven with Florschuetz's own flow (row 1,
    Re_j = 1e4) through a `ChannelElement`:

    | | bias vs Fig. 5 | MAE |
    |---|---|---|
    | before | +268% (range +45% to +762%) | |
    | after | -0.2% | 5.2% |

    The correlation alone gets 5.6% MAE. Fig. 5 (27 points, absolute
    first-row Nu1) is newly digitised; Fig. 6 is Nu/Nu1, where Re_j
    cancels.
  - **Elements now pass their flow area** to `ConvectiveSurface.htc_and_T`
    (new `flow_area` argument), so a non-circular channel with a separate
    `Dh` recovers its mass flow exactly. The fallback
    `pi Dh^2/4` is off by `(Dh/diameter)^2`, -79% on the validation
    channel.
  - **Results change wherever the convective area differs from the
    channel cross-section.** GUI defaults (a 1 m x 0.05 m channel), single
    jet at 0.01 kg/s: Re 1.4e7 -> 1.7e5, h 84,584 -> 3,027 W/(m^2 K).

- **Wall-coupling `dh_dmdot` was wrong for ribbed and impingement surfaces
  (#456).** It was off by a factor of 3.3 (ribbed) and 20 (impingement) in
  the cases checked. `dh_dmdot` is a derivative with respect to the
  element's own mass flow (what the solver's wall coupling multiplies it
  by), but:
  - **Wrong area:** the ribbed and both impingement paths took the mass flow
    on the CONVECTIVE area instead of the channel cross-section, scaling the
    derivative by `A_cross/A_surface`.
  - **Missing term:** the ribbed side `h_r = St rho v cp` differentiated only
    `St`, dropping the dominant `St d(rho v cp)/dmdot` term.

  | case | before | now = finite difference on the channel mass flow |
  |---|---|---|
  | ribbed (the issue's case) | 2320 | 7652 |
  | impingement array (surface area 0.05 m^2) | 3466 | 353,081 |

  - **Jacobian only:** converged solutions are unchanged; a pressure-driven
    coupled ribbed network reaches the same temperatures within 3e-6 K.
  - **New FD tests:** ribbed (1/2/4 walls, two areas), impingement array and
    single jet. The old impingement tests converted through the convective
    area, so they agreed with the bug.

- **Channel heat transfer and friction no longer jump at regime switches
  (#448).** `channel_smooth`, `nusselt_circular_channel` and
  `htc_circular_channel` (gnielinski branch) switched regime with hard `if`s:
  Nu jumped +80% and f +64% at Re 2300, and rough-wall friction stepped from
  smooth Petukhov to Colebrook at Re 4000 (+42.7% at e/D 0.02). Each switch
  is now a C1 smoothstep, and each regime stays EXACT outside its band:

  | band | blend |
  |---|---|
  | Re 2300-3000 | f `64/Re` -> turbulent; Nu laminar (3.66/4.36) -> Gnielinski |
  | Re 3000-4000 | rough walls: clamped Petukhov -> Colebrook |

  Laminar flow keeps its exact Nu and `64/Re` up to 2300 -- the band starts
  where transition does, and ends where Petukhov's stated validity begins.
  **What moved:** only Re inside the bands (outside them, at machine
  precision: |change| < 1e-15). Inside, Nu and f sit between the two
  regimes they join -- e.g. Re 2450: Nu 7.35 -> 4.10, f 0.0456 -> 0.0284
  (the old values were Gnielinski and Petukhov already at 2301). The
  validation scorecard is byte-identical.

  **The Jacobian is exact where it was not.** `channel_smooth`'s Gnielinski
  `dNu/dRe` and Colebrook `df/dRe` were central differences; its Petukhov
  `df/dRe` ignored #449's clamp, and the Gnielinski and Petukhov stencils
  used smooth Petukhov even for rough walls. All now come from the same
  functions as the values: `friction_channel_and_derivative`,
  `nusselt_channel_gnielinski_and_derivative`, and an implicit
  `friction_colebrook_dRe`. `channel_smooth`'s units entry named a `T_wall`
  argument that does not exist and omitted six others; corrected.

- **Smooth-pipe Gnielinski is C1 at Re 3000, and its Jacobian is exact
  (#446).** Every smooth-pipe path held Petukhov friction with
  `max(Re, 3000)`, which left `dNu/dRe` jumping by 25% at Re 3000
  (value continuous) -- a Jacobian kink for any element in the transitional
  range. The hold is now `friction_petukhov_clamped`, a C1 quadratic soft-max
  over Re 2500-3500, used at all five sites (`nusselt_gnielinski`,
  `nusselt_petukhov`, and the three channel paths).

  **What moved:** nothing outside Re 2500-3500 (bit-identical); inside,
  Nu by at most -1.53%, at Re 3000 (friction there -1.36%). The validation
  scorecard is byte-identical before and after.

  `nusselt_and_jacobian_gnielinski` now returns the exact derivative
  instead of a central difference, and for Re <= 0 the laminar value with
  zero slope (flagged extrapolated) instead of `(0, 0)` INVALID.
  `nusselt_gnielinski_with_derivative` / `_smooth_with_derivative` expose
  the analytic form in C++; the ratio-set handover (#444) uses it, replacing
  its private copy.

- **Taslim & Spring's figure 11 is on its printed axis, and its `e/D` are
  numbers.** The axis is `Reynolds No. (x10^-4)`, so the 16 figure 11
  series carry `x_scale: 1.0e4`; 23 `geometry.e_D` entries ended in a stray
  full stop (`0.083.`) and parsed as strings.

- **Figure 4.51 is scored at the rig its paper states.** Han, Zhang and Lee
  (1991) p. 590 gives square, `e/D` 0.0625, `P/e` 10; the 18 measured
  series carried `W/H` 1 or nothing, so the runner probed
  `han_park_1988_angled` at its validity midpoint `P/e` 15. Eq. 4.18's
  square-channel `G` carries `(P/e/10)^0.1`, so that read the model ~4%
  HIGH -- flattering it. In-domain bias moves 60 deg parallel -6.8% ->
  -9.3%, 45 deg -13.8% -> -17.2%; the 90 deg classes do not move
  (`han_1988_orthogonal`'s `G` has no geometry term). Both parallel classes
  now read low by more, the same known Eq. 4.18 accuracy limit. The figure
  4.51 test keeps its bounds on the in-domain cells and holds only the sign
  on the single extrapolated points, which the corrected geometry pushed
  just past them (60 deg 15.2%); `test_figure_451_carries_its_stated_rig`
  pins the geometry, since a regression here would IMPROVE the score.

- **Figure 4.51's disputed 45 deg classes resolved against the primary
  paper (#403 item 1).** Han, Zhang and Lee (1991) JHT 113, 590 -- now on
  disk -- prints `G` and `G_bar` fits for all nine of its rib
  configurations (Table 2). Every confirmed figure 4.51 class lands on its
  own panel's fit within 2.8%; for 45 deg parallel and crossed, BOTH panels
  land on the printed `G` (bias +3.8% / -2.1%) and miss `G_bar` by 13-19%.
  The `G_bar` panel carried the `G` marks twice, so `G_45par`/`G_45crs` are
  now confirmed and `Gbar_45par`/`Gbar_45crs` stay disputed and unscored.
  The "impossible `G_bar/G` < 1" was the duplicate, not physics: the
  printed ratio is 1.07-1.26 for all nine shapes.

  **What it surfaced.** 45 deg parallel `G` now scores against
  `han_park_1988_angled` and reads -13.8%: the recorded "Eq. 4.18 reads
  low" accuracy limit, now in a second independent dataset. The test that
  asserted every figure 4.51 parallel class within 15% is rewritten to pin
  that as the accuracy result it is (bias must stay negative), not
  loosened to pass.

- **Figure 4.53's smooth-friction reference is the one the source prints.**
  The same paper states `f_0 = 0.046 Re^-0.2` (its Eq. 2) for the rig whose
  continuous-rib data figure 4.53 reprints; #440 had carried CR-4015's
  four-sided `0.079 Re^-0.25` across before the paper was on disk. MAE
  7.9% -> 7.6%. It also states the rig (no longer borrowed from 91-GT-3),
  and that its Nu ratio is averaged over X/D 0-20 INCLUDING the entrance
  region -- now declared on the series and in `conventions.py`, since
  `han_1988_orthogonal` is fully developed.

- **CI skipped every test job for changes under `validation/`, and the
  status check passed anyway.** The `source` path filter listed
  `src/ include/ python/ tests/ scripts/ examples/` but not `validation/`,
  so a pull request touching only that tree ran the two path detectors and
  the two aggregators and nothing else. PR #437 merged with **10 of 14
  checks skipped** that way.

  Two kinds of change were affected, and the second is worse:

  - **the scored dataset** -- `validation/cooling/data/**`, which
    `test_every_digitised_series_passes_its_figure_card` and the scorecard
    tests read directly. Those are precisely the tests a bad dataset
    breaks.
  - **the runners** -- `scorecard.py`, `conventions.py`, `fidelity.py` and
    the per-physics runners live under `validation/`, not `python/`. A
    change to any of them skipped CI unless it happened to touch
    `python/tests/` as well, which is the only reason #432 and #433 were
    exercised.

  `validation/**` added to the `source` filter in both `ci.yml` and
  `validation-tests.yml`. A dataset-only change now also triggers the C++
  matrix, which is wasted but harmless; splitting the filter further would
  not have covered the runner case.

### Security

- **combaero-gui declares floors for its web stack.** Before this, `fastapi`, `uvicorn`, `pydantic` and `pandas` had no lower bound. An environment that already had an older FastAPI or Starlette kept it, including versions with published advisories. The new floors clear every known advisory (OSV, checked 2026-10-25):
  - `starlette>=1.3.1`: Host and path poisoning, form-limit bypass, HTTP method dispatch, and StaticFiles SSRF on Windows (CVE-2026-48710/48817/48818/54282/54283). Starlette is a transitive dependency, declared directly so that the floor binds.
  - `h11>=0.16.0`: malformed chunked bodies (CVE-2025-43859).
  - `fastapi>=0.133.0`: the first release that allows Starlette 1.x.
  - `pydantic>=2.7.0`, `uvicorn>=0.30.0`, `pandas>=2.2.3`.
  - The development lock moves Starlette from 0.49.3, which was affected, to 1.10.0.

- GUI CSV export: a text cell (a label, an id) that starts with `=`, `+`, `-`, `@`, a tab or a carriage return gets a leading `'`, so a spreadsheet opening the export shows it as text instead of running it as a formula. Labels travel with shared network files; `=HYPERLINK(...)` used to be written as is.

### Documentation

- **The effusion gas-side augmentation Andrews 88-GT-290 requires is tabulated against both baselines,** for both plates, with and without the film (`validation/cooling/extractions/andrews_effusion_overall_eta.md`). The ranges quoted in the `EffusionPlateElement` docstring and the API reference were stale against current main; they are corrected and now point at the table.

- **Florschuetz Eq. 8's pressure range is documented.** `crossflow_to_jet_ratio_at_row` is incompressible, and its measured range reaches a channel pressure of 0.72 of the jet plenum's (NASA CR-3630 Table D.1). There it matches the measured distributions to 3.4% RMS, the same as the compressible network chain. Below about 0.7 it is extrapolation.

- **`BorderCarnotLossElement` is not for junction ports (#272).** The API_PYTHON example put one on a `MultiPortChamberElement` lateral port. The Mynard closure already carries that turning loss, so the pair double-counts it.
  - **Measured on Bassett's 90 deg, psi = 1 dividing tee:** K6 = 0.867 from the closure alone (Bassett 0.867), but 1.249 with the element at q = 0.5. That is +16% at q = 0.3 and +78% at q = 0.7.
  - **Docs corrected:** the example, the docstrings and the C++ header comments.
  - **Validation status:** the element is described as never validated. Its only check was paired with the junction model removed in 0.6.0.

- **Junction Tier 2 validation (#272).**
  - **Wang 2014** is scored as cross-source accuracy: 116 of 200 points, MAE 0.108, with error flat in Mach. It is the only junction source that set no constant in the closure.
  - **Network scorecard:** prints each source's role (`selected`, `tuned`, `x-source`).
  - **Perez-Garcia 2010 is declared not a validation source.** Its K_hat depends only on the branch Mach numbers, and a junction with no loss reproduces its Table 1 within U95 in 57 of 72 cells.
  - **Its Table 1 transcription is corrected:** seven values were wrong, as were the K_hat definition and the q axis in the docs.
  - **Wang's "4% uncertainty band" is removed:** the paper states none.

- **The Han rib sets are scored on an independent lab:** NASA CR-4396 Vol. II (Johnson, Wagner & Steuber 1993, Pratt & Whitney), the first rib Nusselt data from outside Texas A&M, labelled cross-source accuracy. 48 digitised points from Figs. 6.1/6.2, first leg (`validation/cooling/data/johnson1993/`).
  - On a common Dittus-Boelter basis `han_1988_orthogonal` reads the normal trips +12.1% high and `han_park_1988_angled` the 45 deg skewed trips +18.4% high: a flat level, no Re trend. Every point is above both sets' e/D range.
  - The source's Kays-Perkins reference goes through a conversion registered in `validation/cooling/conventions.py` (0.873 at Pr 0.72). Compared raw, the same sets would read within 3%.

- **Why the junction energy checks stay next to the second-law check (#598).** `MultiPortChamberElement`'s docstring now says why its energy check is not replaced by sigma. At a junction whose ports differ in temperature, sigma's thermal mixing term can mask a mechanical-energy gain that no passive junction can make. Where the ports are isothermal, the two checks agreed on every audited root.

- Ejector validation scores the Huang kernel on the Florence R245fa chiller's third duct (Mazzelli & Milazzo 2015, Tables 3-4; A_3/A_t 9.72 on the same nozzle as Eames 2013). It is cross-source, with tabulated measured mass flows giving the first entrainment ratio with a measurement band on this rig. The kernel's omega is 25-30 % low (Eames 2013's COP-derived omega had it about 2x high); P_c* is bracketed by the two closures (Huang -18 %, Kracik & Dvorak +13 %); a pre-registered interpolation of #687's duct trend was missed and is reported (#697).

- Junction docs explain why a merging tee with equal supplies and equal upstream legs splits exactly 50/50 at any angle: every joining closure gives the suppliers one junction static pressure, which Blaisdell & Manson's data support to within their scatter. They also warn that `ConstantKTeeElement.bassett_design_K` can return a net-negative joining pair at shallow angles, which the second-law check then rejects (#779).

- **Acoustics API examples that did not run.** `docs/API_PYTHON.md` called `tube_axial_modes(..., bc1=, bc2=)` and `orifice_impedance_with_flow(mdot=, area=, ...)`, neither of which exists; `docs/API_CPP.md` gave `tube_axial_modes` a `std::vector<double>` return. The `can_annular_eigenmodes` docstring pointed Python users at `can_annular_eigenmodes_with_method`, which is not bound.

- **The combustion API references match the code.** `docs/API_PYTHON.md` showed `cb.complete_combustion(state, smooth=True)`, which raises `TypeError`. The real call is `(T, X, P, smooth_phi0, smooth_phi1, k0, k1)`. `docs/API_CPP.md` listed the old `bool smooth` signatures and a `PressureLossContext` with fields it does not have. Both now give the real signatures, the rich-mixture model (all fuels burn 1/phi, no CO/H2) and what the phi = 0 smoothing costs: n_req = a - ln2/k0 at every lean phi, about -0.3 to -0.6 K on T_ad.
  - `units_data.h` labelled `oxygen_required_per_kg_*` and `dryair_required_per_kg_*` as mol/kg. They return kg/kg (CH4: 3.99 kg O2/kg).

- Combustor docs brought in line with the code (found by a combustor audit):
  - The "Network Nodes" example now runs. It had unknown constructor keywords on every node, and a `MomentumBoundary` that does not exist.
  - The (f, J) interface example now names the real `combaero._solver_tools` signatures.
  - API_CPP gives the combustor functions' real `MixerResult` signatures.
  - `EnergyBoundary(fraction=...)` is documented as a fraction of the sensible enthalpy, and of the products' on a combustor (#481), as the code has done since.
  - The combustion smoothing comments state the actual constants.
  - A failed equivalence-ratio diagnostic now raises a `RuntimeWarning` instead of printing `DEBUG:` to stderr.

- **The ejector kernel is scored on a refrigerant rig:** Eames et al. 2013, an R245fa chiller with two mixing ducts on one nozzle (`validation/ejector/eames2013_check.py`).
  - Measured critical pressure: Huang's printed chain (the default) is 14% low and Kracik & Dvorak 20% high, and the two bracket every point.
  - The larger duct moves both closures by 10-22 points. The miss is geometry in the shared front end, not the mixing closure.
  - The entrainment ratio is scored against the COP the authors themselves doubt, with that doubt recorded.

- **The ejector's operating anchors are scored cross-source** (`validation/ejector/anchor_scorecard.py`), on Henry 2007's measured air ejector and Akbarnejad & Ziabasharhagh 2025's steam CFD. The dead-head back pressure P_b0, the far end of the subcritical droop, runs 78-274% above every rig's breakdown, and the error grows with A_3/A_t. That is why deep-subcritical roots at large mixing-area ratios are refused by the second law (#598). The README now declares the validated envelope of each closure. No model change.

- Ejector validation: Chen et al. 2013 (Figs. 5a and 7a) digitised. It adds a third measured rig (a large rectangular air ejector with error bars) to the cross-source anchor scorecard. Chen's own model curves are a second author's evaluation of Huang's chain. On the new rig the candidate sudden-expansion dead-head lands 3-5% above the breakdown, against 62-68% for today's anchor.
- The scorecard now also scores a candidate pair (Huang's printed P_c* chain with the sudden-expansion P_b0) over whole measured curves: 12-23% RMS against 29-59% for today's closures. The pair is not adopted yet.
- `validation/ejector/README.md` documents why the printed chain misses both authors' own numbers. Fourteen single-change hypotheses are tested, along with unit slips, real-fluid R141b and property-state sets, and none explains the gap. Henry 2007's primary pressures are shown to be gauge, and Chen's rectangular rig's absolute.

- The ejector validation notes cover three more sources:
  - Bouhanguel et al. 2011: an air ejector's observed flow regimes, which only the new default P_c* matches without an edge-case reading (`validation/ejector/bouhanguel2011_regime_check.py`);
  - Lucas & Koehler 2012: transcritical CO2, outside the ideal-gas element, whose choke is 34-51% below the measured motive flow;
  - Liao 2008: no data with a back-pressure sweep.
- The validated envelope now also lists two limits: a dense or flashing motive fluid, and a mixing tube too short to finish mixing.

- Cooling validation: every scored rib series is now scored by every eligible correlation set, not only CR-4396. Eligibility tightened: another set qualifies only on a declared or registered-conversion quantity, and a single-angle set never scores a series that does not state its angle. 47 series gain rows. The Han, Zhang and Lee (1991) sets gain fidelity rows on their own reprinted Fig. 4.51 data (+0.4% to +5.5% bias). `han_park_1988_angled` and `han_1989_narrow_channel` share a citation but are fidelity only on their own figure (4.47 and 4.48).

- **The junction scorecard no longer invents an uncertainty band.**
  - **Before:** `pct_within_uncertainty` applied a flat 0.05 to every source. Bassett's metadata said 0.05, Hager's 0.05 and Idelchik's 0.10, none with a source.
  - **The full texts state no measurement uncertainty.** Bassett and Hager compile other experimenters' data, and Idelchik's preface disclaims accuracy outright.
  - **Now:** each paper's band is null with its reason, and "within" reads `-` rather than a score.
  - **Idelchik provenance:** its theta = 30 tables (diagrams 7-1/7-2) are calculating formulas with no experiment behind them, per the handbook's own diagram list. This is recorded in the metadata and README.

- Cooling validation: a series can be scored by several correlation sets (`also_scores`), each on its own scorecard row. The list must equal every set the harness can score the series with (`runner.eligible_sets`, by rule), so no set can be left off because it reads badly. The report adds a "scored by several sets" section that lists sets in name order, ranks none, shows the spread of bias across them, and flags when all of them come from one lab. NASA CR-4396 now scores all six eligible rib sets: Fig. 6.1 reads +5.2% to +12.1% bias, Fig. 6.2 -1.6% to +18.4%, and Rallabandi 2009 reads -0.6% on its 12 in-box points.

- Documented the 1000 K species-polynomial join in the C++ API reference: values are continuous, slopes one-sided (dk/dT -0.01% air, -0.1% lean CH4 products). The data stay unblended, and a test pins the bound.

- The combustor's `"equilibrium"` method is documented as what it is: a PARTIAL equilibrium (complete combustion, then reforming and water-gas shift among CO, H2, CO2, H2O), with no dissociation by design. It equals `"complete"` when lean and is cooler when rich. The earlier "physically correct for rich or high-T conditions" claim is removed: without dissociation it is not the high-T equilibrium. The decision is recorded at `CombustionMethod` in `include/state.h`.
- GUI: the method selector reads "Complete (CO2 + H2O)" and "Partial equilibrium (CO/H2 shift)", with a line of help text under it. The stored value stays `"equilibrium"`, so saved networks are unchanged.

- `pressure_loss.py` no longer claims its correlations give "exact analytical Jacobians" through the C++ combustor chain. The loss is evaluated on a `PressureLossElement` in Python: the theta column is analytic, the other dependences are central differences, and the C++ chain receives a zero loss. That is by design for now (#718). The callable is a prototype interface for user-supplied loss maps (theta and fuel splits) whose input design is still open.

- A source register, `validation/sources.yaml`, with one entry per local reference document (126): citation, content, what it was read for, the verdict, and the issues, records and code that use it. `scripts/sources.py` searches it, lists by status or topic, and checks it (schema; unregistered local documents).

- Cooling validation: Taslim & Spring (1987) two-wall rib friction at AR 0.5 and 3.5 is now scored against `han_1988_orthogonal`, not only the square channel. An independent lab's data, outside Han's range, labelled as accuracy: AR 0.5 reads -12% to +13%, AR 3.5 (Han W/H 0.286) about 38% low. A test now fails if any measured rib series that a set can score is left unscored without a recorded reason.

## [0.6.0] - 2026-09-09

### Changed
- **The momentum-CV junction drops its version suffix.** `MPCEv2Element` is now
  `MultiPortChamberElement`, taking the name freed by the removal below. The
  abstract base that owns the port machinery is now `MultiPortChamberBase`; it
  has three subclasses -- the junction, `ConstantKTeeElement` and
  `EjectorElement` -- so it stays a class of its own rather than being merged
  away. `MultiPortChamberElement` and `ConstantKTeeElement` are now exported
  from `combaero.network`, which they were not before.

  Renamed with them: the module `mpce_v2_element` -> `mpce_element`, the
  offline Jacobian cross-check `_mpce_v2_jacobian` -> `_mpce_jacobian`, the C++
  entry point `mpce_v2_residuals_and_jacobian` ->
  `mpce_residuals_and_jacobian`, the validation adapter `MPCEv2Network` ->
  `MPCENetwork`, and the `test_mpce_v2_*` files.

  **Breaking for anyone importing the old names.** `from combaero.network
  import MultiPortChamberElement` still works and now yields the concrete
  junction rather than an instantiable base; code importing `MPCEv2Element`
  from `combaero.network.mpce_v2_element` must move to
  `combaero.network.mpce_element`. The GUI's persisted `mpce_tee` node type is
  UNCHANGED, so saved networks still load. No behaviour change: the junction
  scorecard is identical at 2178 of 2546 converged and a mean absolute error of
  0.4809.

  Past CHANGELOG entries and `docs/archive/` keep the old names, because they
  record what the code was called at the time.

### Removed
- **`MultiPortChamberBase`'s own junction model, and the C++ stack behind
  it.** Deprecated in 0.5.0, removed here as announced. The class remains as
  the ABSTRACT BASE owning the topology and port machinery that
  `MultiPortChamberElement` and `ConstantKTeeElement` inherit -- port ordering, the
  inlet/outlet sign map, area and angle resolution, the wiring checks, and
  `diagnostics`, which reports the port map and computes no physics. What went
  is `residuals` and `verify_solution_consistent`, so the class is no longer
  instantiable on its own.

  Also removed: `multi_port_chamber_residuals_and_jacobian` (C++ and its
  pybind11 binding), the `MultiPortChamberResult` and `PortImpulseJacobian`
  structs, the `combaero._solver_tools` re-export, and the `units_data.h`
  entry. That function had exactly one caller -- the method removed above.
  `include/multi_port_chamber.h` is KEPT: it also holds `border_carnot_L` and
  `HAGER_FRACTION`, which `border_carnot_loss_residual_and_jacobian` is built
  on, and `BorderCarnotLossElement` is unaffected.

  Migration: use `MultiPortChamberElement` (the Mynard closure) or `ConstantKTeeElement`
  (fixed handbook K). On the same Bassett separating cells the removed model
  scored a mean absolute error of 0.5260 against `MultiPortChamberElement`'s 0.0564 and
  converged on 77 of 105 points against 94; it was also energetically
  inconsistent for joining flow and sign-symmetric, so it admitted mirror
  roots. Nothing in the package or the GUI instantiated it. Pin
  `combaero~=0.5` to keep running a network built on it.


## [0.5.0] - 2026-09-09

### Added
- **`NetworkSolver.solve` now reports why a solve ended, separately from
  whether it succeeded.** Five new keys: `__converged__` (the root finder's own
  verdict), `__consistent__` (the elements' physical-consistency verdict, with
  **`None` meaning not checked** -- never `True`), `__inconsistent_elements__`,
  `__outcome__` (a `SolveOutcome` `StrEnum`, always set) and
  `__worst_residuals__` (the rows carrying the residual, largest first).
  `__success__` keeps its exact meaning, converged **and** consistent, so
  existing callers are unaffected. The motivation: `__success__` is False both
  when Newton never got there and when it found a root a junction then rejected
  as unphysical -- the latter reporting a tiny `__final_norm__` next to
  `success=False` -- and everything else about the failure was only available
  by matching substrings of `__message__`, part of which comes from SciPy.
  `SolveOutcome` is exported from `combaero.network`. The GUI's `NetworkResult`
  surfaces the same fields, and its docstring no longer describes `success` as
  "converged within the residual tolerance", which it has not been since the
  consistency checks landed.
- **The tuned constants in the MPCEv2 junction closure are now named,
  labelled with their provenance, and switchable for measurement (issue
  #271).** `combaero.network._mynard2010` gains `MYNARD_ETA_A0 = 0.8`,
  `MYNARD_ETA_A1 = -0.2` (Mynard & Valen-Sendstad 2015 Eq 36, fitted to their
  Fig 4 CFD at blood-flow Reynolds numbers) and `FLOW_RATIO_DAMPING = 0.02`
  (a regulariser from the Matlab reference, not in the paper), each documented
  at its definition. `MultiPortChamberElement` and `junction_loss_coefficient` accept a
  new `eta_scale` keyword (default `1.0`, the faithful port; `0.0` switches
  Mynard's energy-transfer term off) so the term can be scored on and off
  against the digitised validation data alongside `joining_etransfer_alpha`.
  **No behaviour change at the defaults** -- the closure is bit-identical.
  The provenance label on `DEFAULT_JOINING_ETRANSFER_ALPHA` is rewritten: the
  value was fitted with a mirrored Bassett K11 anchor (#283) and the
  "145-point independent validation" it cited is not in the script that
  produced it; `scripts/calibrate_mpce_v2_etransfer.py` (stale import path,
  absolute data path, no validation code) is retired in favour of
  `validation/junction/calibrate_etransfer.py`.
- **Ejector operating-regime extension wired into `EjectorElement` (Phase 2).**
  The network `EjectorElement` now resolves all three physical regimes in one
  C1 Newton-solvable residual system instead of the critical (double-choked)
  plateau only: critical, subcritical droop, and the unchoked-primary subsonic
  jet pump. A previously non-converging network -- a low primary mass flow
  (below the choke threshold at its back pressure) drawing suction against
  atmospheric -- now converges to the physical jet-pump root (primary Pt floats
  to ~the back pressure, `dp >= 0`, entrainment ratio `omega ~ 7`) instead of
  failing or landing on a spurious double-choked root with the primary pressure
  below the back pressure. The regime is selected smoothly by two smootherstep
  weights (`s_choke` on `mp / mdot_choked`, `s_sub` on `outlet.Pt / P_c*`) with
  no branching discontinuity, and the whole 4-row residual system plus its
  analytic Jacobian are assembled in C++ (a single
  `ejector_element_residuals_and_jacobian` chaining the scalar closures through
  a forward-mode `DualN<9>`), matching the whole-element `(f, J)` practice of
  the base `MultiPortChamberBase`/`TeeJunctionElement`; the Python element is
  a thin relabeling shim. The critical regime is exactly the limit of the new
  blend so existing critical-mode behaviour is unchanged.
  `diagnostics()` now reports the ACTUAL entrainment ratio (`omega = ms/mp`,
  regime-independent), a `critical_mode` regime flag, the mixing-plane pressure
  `p_py_pa`, and suppresses the critical-only `p_c_star_pa`/`omega_critical`
  when the primary is unchoked. `verify_solution_consistent` no longer demotes
  subcritical/jet-pump solutions (all regimes are now modeled). The GUI ejector
  warm-start seeds the jet-pump operating point when the primary cannot choke.
  See `validation/ejector/OPERATING_REGIMES_DESIGN.md` sec 6c/8.3.
- **Ejector unchoked-primary jet-pump reference closure (Phase 1B of the
  operating-regime extension).** `jet_pump_entrainment_ratio` in
  `combaero.network._ejector_huang1999` (re-exported from
  `validation/ejector/models/huang1999.py`) models the subsonic jet-pump regime
  the double-choked construction cannot: a Keenan-Neumann-Lustwerk constant-area
  mixing solve (both streams expand to a common mixing pressure; the physical
  mixing pressure recovers to the back pressure). The mass+momentum+energy
  mixing reuses the existing Kracik & Dvorak closure with a subsonic primary
  core (verified bit-identical to a direct momentum balance), so the jet-pump
  and ejector branches share one mixing model. Validated to reduce to the
  lossless incompressible jet-pump relation N(M,R) (Sanger, NASA TN D-4445) as
  Mach -> 0 (error -> 0.03%), and to give a bounded omega O(few) on the reported
  degenerate network instead of the ~33 the choked construction returns.
  `recovery_efficiency` stays 1.0 (unfitted) but is audited against Sanger's
  measured data as a ~1.4x overprediction (realistic ~0.7). The golden fixture
  gains an air jet-pump back-pressure sweep (omega + d(omega)/d(P_b)).
  Reference-model only -- not yet wired into `EjectorElement` (Phase 2). See
  `validation/ejector/OPERATING_REGIMES_DESIGN.md` sec 8.2 for the provenance
  record.
- **Ejector subcritical-droop reference closures (Phase 1A of the
  operating-regime extension).** Three closed forms in
  `combaero.network._ejector_huang1999` (re-exported from
  `validation/ejector/models/huang1999.py`): `dead_head_back_pressure` (the
  zero-entrainment upper anchor P_b0 -- Kracik & Dvorak's *existing* mixing
  closure at omega = 0, no new calibration constant),
  `subcritical_entrainment_ratio` (Tier-1 linear droop between P_c* and P_b0),
  and `blended_entrainment_ratio` (C1 smooth-min joining the critical plateau
  to the droop through P_c*). `critical_back_pressure` was refactored onto a
  shared y-y-state + mixing helper so its 31 golden values regenerate
  byte-identically. The Huang golden fixture now also carries per-case P_b0 and
  a subcritical back-pressure sweep (omega_eff + d(omega_eff)/d(P_b)) as the
  Phase 2/3 target. Reference-model only -- not yet wired into `EjectorElement`
  (that is Phase 2). The near-linear droop is validated against measured air
  data (Henry et al. 2007 HEFAT2007 Fig. 5, R^2 >= 0.94); see
  `validation/ejector/OPERATING_REGIMES_DESIGN.md` sec 8.1 for the provenance
  record.
- **`EjectorElement`: production supersonic-ejector network element**
  (`combaero.network.ejector_element`), critical (double-choked) mode only.
  Reuses Huang (1999)'s closed-form entrainment ratio (Eqs. 1-8) together
  with Kracik & Dvorak (2016)'s lambda/q(lambda)/z(lambda) mixing closure
  for the critical back pressure (their Eqs. 7-13), replacing Huang's own
  original momentum+shock+diffuser P_c* chain (Eqs. 9-18) in the reference
  model -- validated against Huang's 31-row R141b dataset (mean/max error
  6.2%/13.2%, vs. 25%/35% for Huang's original chain and 9.0%/16.0% for a
  Lienhard/McGovern zero-loss closure also evaluated; see
  `validation/ejector/README.md` for the full comparison and paper
  references). Sits on the `MultiPortChamberBase` 3-port junction
  topology (primary/secondary inlets, one outlet); `P_jct` is repurposed
  as a diagnostic critical back pressure P_c*, not the outlet's actual
  pressure, since critical-mode entrainment ratio is independent of back
  pressure by definition. Supplies its own port-face areas from geometry
  (primary = nozzle exit, secondary = mixing - nozzle exit, outlet =
  mixing), so it resolves even when fed by a bare boundary/connection with
  no channel to infer an area from (unlike a generic tee, which must infer
  or be given `port_areas`). gamma/R evaluated live from combaero's own
  real-fluid EOS (combustion/air-relevant species only) at the
  entrained-flow choking plane. `recovery_efficiency` (default 1.0, no
  artificial loss) is intentionally not fitted to the validation dataset,
  per this project's calibration policy of anchoring on analytical models.
  Residuals and their full analytic Jacobian are evaluated through a C++
  `(f, J)` path (`include/ejector.h` / `src/ejector.cpp`, exposed via
  pybind), per CLAUDE.md's Solver `(f, J)` rule. The underlying
  closed-form physics (`entrainment_ratio`, `critical_back_pressure`,
  `choked_mass_flow`, `EjectorGeometry`) now lives in
  `combaero.network._ejector_huang1999`, following the same
  production-tree placement already used for the Mynard tee-loss model,
  since `python/combaero/` must never import from the dev-only
  `validation/` tree; `validation/ejector/models/huang1999.py` re-exports
  it unchanged for existing validation scripts.
- **GUI: `EjectorElement` is wired end-to-end as a draggable "Ejector"
  element** (Wind icon). Sidebar palette entry, a 3-port node card (primary
  + secondary inlets, one outlet), an inspector form for the three areas
  (throat / nozzle-exit / mixing) and `recovery_efficiency`, drop-time
  defaults, and client-side geometry validation. The backend
  `graph_builder` translates an `ejector` node into an `EjectorElement`,
  reusing the tee's named-port resolution and per-port
  `MomentumChamberNode` auto-insertion (generalized from `tee_ids` to a
  shared `junction_ids`); solved diagnostics (entrainment ratio omega,
  critical back pressure P_c*, gamma) surface through the generic
  `ElementResult` path. `graph_builder` also warm-starts each ejector: it
  traces each port to its feeding boundary, runs the ejector's own
  closed-form physics for a consistent operating point, and seeds the port
  MomentumChamberNode pressures and `P_jct` (only where the user provided
  none) so the stiff double-choked solve converges from cold -- which it
  does not otherwise (no solver init strategy cracks it without the seed).

### Changed
- **The momentum-CV junction's residual and Jacobian are now computed in C++.**
  `MultiPortChamberElement` becomes a shim: it applies the guards -- the degenerate-state
  fallbacks and the wrong-direction soft barrier, which are solver policy
  rather than junction physics -- and then calls
  `_core.mpce_residuals_and_jacobian`, which seeds the whole element over
  `(P_i, Pt_i, outer_mdot_i, Pt_jct)` and returns exact partials for every one.
  This replaces a Python closure evaluation, a sympy-lambdified Jacobian block
  for the canonical separating tee, an N+1-call finite-difference fallback for
  every other topology, and a hand-derived `dR/dP` column.
  **The Jacobian is now more complete, not merely faster.** The hand-derived
  column covered the common port only, but `K` depends on every port's velocity
  and every velocity on its own density; measured against central differences of
  the assembled network Jacobian, the worst relative error falls from 4.4e-3 to
  1.4e-6. Accuracy is unchanged to four decimals on every validation source
  (Bassett 0.1538, Hager 0.0787, Idelchik 0.8880, Wang 0.1184) and convergence
  is slightly better: the scorecard goes from 2108 to 2109 of 2546, and the
  random boundary sweep from 92.9% to 93.7% of draws that admit a root.
  `MultiPortChamberElement.jacobian_method` is retired and now selects nothing; it is kept
  so existing callers do not break. The sympy derivation remains in
  `_mpce_jacobian.py` as an offline cross-check rather than a runtime path.
- **The MPCE junction closure now carries the dividing-streamline pressure
  recovery, and Mynard's fitted energy-transfer factor is off by default.**
  With the transfer factor disabled, the closure's straight-leg coefficient was
  exactly `q^2` -- the plain velocity-difference loss -- at every area ratio and
  branch angle, while Bassett 2001 (K5) and Hager 1984 (xi_t) independently give
  `q^2 - 0.5 q`. The missing `-0.5 q` is the dividing-streamline pressure
  `p* = p + (1/4) rho u^2` acting over the diverted flow fraction; Mynard
  carries the same term but only through the contraction analysis of a turning
  collector, whose control volume degenerates for a collinear one.
  `_mynard2010.DIVIDING_STREAMLINE_RECOVERY` restores it for the single-supplier
  diverging case that both papers analysed.
  Restoring it made Mynard's `eta` a duplicate of the same physics, so
  `MultiPortChamberElement.DEFAULT_ETA_SCALE` is now `0.0`; `eta_scale=1.0` still
  reproduces the faithful port. Measured at pinned operating points on the
  digitised data, the straight-leg RMSE improves from 0.333 to 0.098 while the
  lateral regresses from 0.079 to 0.097: the fitted factor was buying lateral
  accuracy by making the junction a net source of flow work below a lateral
  flow fraction of about 0.25. The closure is now dissipative everywhere and
  reproduces the analytical identity `K_lateral - K_straight =
  q (0.5 - 2 cos(0.75 theta)) + 1` at equal areas. The analytical Jacobian was
  re-derived with the term and with `eta_scale` as a parameter rather than
  baked in at 1. **Junction pressure drops will change**, most in dividing flow
  at a small lateral fraction.
- **Flow-area inference now consults every area-bearing neighbour, and refuses
  to invent an area it cannot determine (issue #262).** Seven `resolve_topology`
  sites searched only for a neighbouring `ChannelElement` and, failing that,
  substituted a hard-coded constant (0.01 / 0.02 / 0.1 m^2). With a non-channel
  neighbour -- a momentum chamber, a combustor, an ejector outlet -- that
  fabricated a geometry the network does not contain: the reported case
  (`momentum_chamber (0.15 m^2) -> area_change -> channel (0.636 m^2)`) resolved
  `F0 = 0.01`, a 63.6x expansion out of nothing, and stalled the solver with a
  message that blamed the solver rather than the geometry. In simpler chains it
  is worse than a stall -- it converges, inflating the chamber stagnation
  pressure by 4% at 1 kg/s, 64% at 5 kg/s and 342% at 20 kg/s.
  All seven sites now share one inference path (channels first, so networks that
  resolved before resolve identically; then other area-bearing elements; then
  adjacent nodes whose area was set explicitly). An auto-sized area is never
  used as a source.
  **Breaking:** where the area *is* the geometry being modelled --
  `AreaChangeElement` (`F0`/`F1`), `ChannelElement` (`diameter`),
  `TeeJunctionElement` (`F_C`), `BorderCarnotLossElement` (`area`) -- an
  unresolvable area now raises a `ValueError` naming the component, the
  parameter that clears it, and where the search looked, instead of silently
  defaulting. Networks that relied on those defaults must set the value
  explicitly or give a neighbour a known area; both clear the error. Where the
  area is only a dynamic-head reference (`MomentumChamberNode`,
  `CombustorNode`, `PressureLossElement`) a nominal default still stands, since
  a network can legitimately never read it -- `PressureLossElement` warns when a
  head-loss correlation will actually consume the defaulted value, and stays
  quiet for fraction-based correlations that never touch it.
- **`EjectorElement` now evaluates its residuals and a FULL analytic
  Jacobian (all four rows) through the C++ `(f, J)` path** rather than the
  pure-Python reference physics plus a finite-difference Jacobian fallback.
  The ejector closed forms and their analytic derivatives (added to
  `include/ejector.h` / `src/ejector.cpp`) are exposed via pybind
  (`combaero._solver_tools.ejector_*_and_jacobian`). The element's Jacobian
  is exact w.r.t. the four Newton-solved thermodynamic inputs (primary and
  secondary stagnation pressure and temperature); `gamma` and `r_gas` are
  frozen coefficients (recomputed each residual call, so the converged root
  is exact -- only the weak implicit `d(gamma)/d(Tt)` term is dropped from
  the step). `combaero.network._ejector_huang1999` is now the validation
  reference only.
- **`MultiPortChamberBase` subclasses can now override
  `row_scale_kinds()`** to declare their own per-row scale-kind pattern
  ("p" or "mdot") for `solver.py`'s row-scaling vector, instead of the
  solver unconditionally assuming every subclass shares the base class's
  N-pressure-rows-then-one-mass-row pattern. Needed for `EjectorElement`,
  whose four rows are the opposite pattern (three mass-flow rows, one
  pressure row); the default matches prior hard-coded behaviour so
  existing subclasses (`ConstantKTeeElement`, `MultiPortChamberElement`) are
  unaffected.
- **Dev tooling: Biome pinned to 2.5.10 across all three sources.** The
  Dependabot frontend bumps of `@biomejs/biome` (2.5.6 -> 2.5.9 -> 2.5.10 in
  `gui/frontend/package.json`) are matched by the `biome.json` `$schema`
  and the `.pre-commit-config.yaml` hook pin, per the frontend
  version-sync checklist.

### Deprecated
- **`MultiPortChamberBase`'s own junction model, for removal in 0.6.0.**
  Its `residuals`, `diagnostics` and `verify_solution_consistent` now emit a
  `DeprecationWarning` once per element. Use `MultiPortChamberElement`: on the same
  Bassett separating cells it scores a mean absolute error of 0.0564 against
  this model's 0.5260 and converges on 94 of 105 points against 77, and it
  computes its residual and Jacobian in C++. This model is also energetically
  inconsistent for joining flow and sign-symmetric, so it admits mirror roots.
  Nothing in the package or the GUI instantiates it.

  **Only the model is deprecated.** The class also owns the topology and port
  machinery that `MultiPortChamberElement` and `ConstantKTeeElement` inherit unchanged --
  13 of its 17 public members -- and that is staying. The warning is guarded on
  the exact type rather than on `isinstance`, because
  `MultiPortChamberElement.diagnostics` delegates here through `super()` and would
  otherwise warn every user of the element that replaces this one.

  0.5.0 is the last release carrying it, so a network built on it can be run by
  pinning `combaero~=0.5`.

### Fixed
- **The junction validation harness built its pressure-driven boundary targets
  from two different operating points.** Bassett Table 1 indexes each loss
  coefficient on the mass-flow fraction in its own leg -- `K5` on
  `mdot_A/mdot_C` (the straight leg), `K6` on `mdot_B/mdot_C` (the lateral) --
  so at one physical state they are `K5(1-q)` and `K6(q)`. All three network
  adapters read both off the same abscissa, which asked the `three_pb` and
  `mfb_two_pb` skeletons for boundary pressures corresponding to a state that
  does not exist. The pairing now lives in one place,
  `bassett2001.separating_pair_at`. Measured: the junction scorecard goes from
  2109 to **2178** of 2546 converged and 2082 to **2152** points scored, with
  the pooled mean absolute error improving from 0.4982 to **0.4809** and
  Bassett's own from 0.1538 to **0.1416**. `mfb_two_pb` changes character --
  it used to drift far from the operating point it was asked for (q=0.8
  settling at 0.857) and now lands within 0.005 across the range, because the
  drift was the solve honestly chasing a mismatched target. Two strict xfails
  recorded as "no root exists, waiting on the K_straight gap" come off: the
  root existed, and the corrected target at q=0.2 is 0.422 rather than 0.122.
  Affects the validation harness only -- no shipped model or solver behaviour
  changes.
- **The junction soft-barrier weight is derived from the network's scales
  instead of being a fixed constant.** `soft_penalty_alpha` multiplies a
  squared mass flow to produce a pressure, so it carries `Pa/(kg/s)^2` and is
  only meaningful against a particular network's scales; the barrier's fixed
  point sits at `slack* = sqrt(dP/alpha)`, and what matters is that as a
  fraction of the flow. Measured on one junction scaled over five decades with
  every dimensionless group held fixed, the weight needed to converge follows
  `1/m_ref^2` exactly -- two decades of weight per decade of size -- and the
  previous fixed value failed on the same junction at a hundredth of its size.
  `NetworkSolver` now derives `alpha = P_ref / (f*m_ref)^2` from the reference
  state it already computes for seeding and hands it to each element that owns
  a barrier, frozen for the solve so the residual and Jacobian are unchanged in
  form. A degenerate reference state falls back to the fixed default rather
  than producing a weaker barrier, and an explicitly set `soft_penalty_alpha`
  always wins. Measured: the junction validation scorecard is identical to the
  digit (2108 of 2546 converged, every source's mean error unchanged) because
  that set already sat on the saturation plateau, while the same junction now
  converges across seven decades of size instead of five.
- **The junction soft barrier no longer has a fixed point that solves park
  in.** `MultiPortChamberElement` replaces the physics with a "soft barrier" when a port
  flows against its declared direction, meant to pull Newton back toward
  `mdot = 0` so a sign flip can restore the declared regime. It could not: the
  penalty is added into the same residual row as the continuity relation,
  `R_i = (Pt_i - Pt_jct) + alpha*max(0, -e_i*mdot_i)^2`, and the element has no
  spare row to give it, so the solver could zero that row by carrying a
  pressure error equal and opposite to the penalty instead of by driving the
  slack to zero. The barrier therefore behaved as a fabricated pressure loss
  with a fixed point at `slack* = sqrt(dP/alpha)`. At the previous
  `soft_penalty_alpha` of 1e7 that put it at 45% of the common mass flow, and
  affected solves plateaued on a residual floor with an ordinary in-regime root
  available. `soft_penalty_alpha` is now 1e11, which moves the fixed point to
  0.45%. Measured: the junction validation scorecard goes from 2080 to 2108
  converged of 2546 with every source's mean error equal or better (Bassett
  0.1548 -> 0.1538, Idelchik 0.8902 -> 0.8880, Hager and Wang unchanged), and
  the random boundary-condition sweep from 98.3% to 99.3% inside the closure's
  documented Mach range. The response is monotone in `alpha` and flat above
  1e9, so this is a saturation point rather than a fitted value. Note that
  `alpha` carries Pa/(kg/s)^2 and is therefore tied to the scales it was
  measured on; a network far outside them wants a scale-aware weight.
- **A failed solve now really does return its best iterate.**
  `NetworkSolver.solve` warns "Returning best iterate" on non-convergence, but
  the returned state and its residual norm were captured immediately after the
  primary root-finding phase. Later phases -- above all the
  Levenberg-Marquardt fallback -- keep evaluating through the same wrapper and
  improve the tracked best, yet only re-pointed the returned state when they
  reached the convergence tolerance, so an improvement that fell short was
  discarded. Measured over 38 non-converged junction solves, 30 returned a
  state worse than the best they had evaluated, by a median factor of 5.8 and
  up to 2e5. The returned state is now re-pointed at the tracked best whenever
  that is closer, before the junction consistency checks run, so those checks
  judge the state actually handed back. **This does not change whether any
  solve converges** -- the junction validation scorecard and the random
  boundary-condition sweep are identical either side of the change.
- **The automatic warm-start retry no longer replaces a closer result with a
  worse one.** When both the primary attempt and the retry failed, the retry's
  result was returned unconditionally. Whichever got closer is now returned,
  with the other's residual norm quoted in the message.
- **The junction validation harness records a real residual norm.** It read
  `__residual_norm__`, a key `NetworkSolver` has never set, so every record's
  residual norm was silently infinity; the solver's key is `__final_norm__`.
- **A network driven only by pressures no longer starts from a hard-coded
  0.1 kg/s.** `NetworkSolver._infer_reference_state` fell back to that constant
  whenever no `MassFlowBoundary` set the scale, so the initial guess was
  independent of the network's size and of the imposed pressure differences.
  The reference flow is now a Bernoulli estimate, `m = A sqrt(2 rho dP)`, using
  the median flow area and the per-element share of the boundary pressure
  spread -- the same estimate `analytical_pt_prop` already made for channels,
  and the same per-element pressure step the guess propagator already used. On
  a sweep of junctions driven by three pressure boundaries, where the implied
  level spans 2e-3 to 36 kg/s and the constant was off by more than a decade in
  27 of 60 cases, convergence rises from 76% to 90%; the junction validation
  scorecard rises from 1697 to 1708 of 2073. Networks that already carry a
  mass-flow boundary are unaffected.
- **The initial guess at a junction now conserves mass and follows the
  boundary pressures.** `NetworkSolver`'s `analytical_pt_prop` seeding applied
  its Bernoulli mass-flow estimate only to `ChannelElement`, and junctions are
  wired with `LosslessConnectionElement`, so every junction port fell back to
  the flat reference flow: on a three-port junction that is 0.1 kg/s in against
  0.2 kg/s out, with a 1:1 split whatever the boundary pressures say. The split
  is now a Bernoulli share of the propagated port pressure differences,
  rescaled so continuity holds exactly at the start with every port in its
  declared direction. The total still comes from the existing topological
  propagator, so an upstream `MassFlowBoundary` continues to set the level and
  no new flow-scale heuristic is introduced. On the junction validation
  scorecard, solves rejected as landing on a physically inadmissible or
  artifact root fall from 227 to 163, and convergence rises from 1625 to 1732
  of 2073, with accuracy unchanged on the records scored under both. The seed
  is opt-in per element class through `seeds_ports_by_pressure_split`, which
  `EjectorElement` sets to false: its port flows follow choking and entrainment
  rather than a pressure split, and it carries its own warm start. Networks
  with more than one admissible operating mode may now report a different one
  than before; see `docs/archive/JUNCTION_OPERATING_POINT_271.md`.
- **`MultiPortChamberElement` now rejects converged states in which the junction would
  create flow work.** Its post-solve `verify_solution_consistent` checked only
  that port flows kept their declared direction; `MultiPortChamberBase` (v1)
  has had an energy check since #229, and v2 had none, so a solution where the
  junction manufactured total pressure was reported as a success. The new
  criterion is a mass-weighted energy balance,
  `sum_in |m_i| Pt_i >= sum_out |m_i| Pt_i`, which is the element's own
  dissipation identity. It deliberately **permits** one branch to gain total
  pressure at another's expense -- that is measured physics in dividing flow,
  and is why Bassett's K5 and Hager's xi_t go negative -- while forbidding a
  net gain, and unlike v1's form it applies to joining flow too. Audited before
  being imposed: the digitised measured curves and Bassett's analytical pair
  satisfy it with margin, while the model violates it below a lateral fraction
  of about 0.2, which is the known `K_straight` gap (#272) surfacing as a
  physical violation. Networks solving a junction in that region will now
  report non-convergence instead of an impossible answer; the solver's message
  names a non-dissipative closure state as one of the causes.
- **`NetworkSolver` results no longer depend on `PYTHONHASHSEED`.**
  `_propagate_pressure_guess` seeded its breadth-first pressure propagation
  from `list(set(...))` over node-ID strings. Set iteration order for strings
  follows `hash(str)`, which Python randomises per process, and every BFS hop
  applies the same pressure decrement -- so the propagated initial guess, and
  with it the Newton basin, changed from run to run on identical input. On
  the junction validation scorecard (issue #271) the `mfb_two_pb` topology
  converged 618 or 629 of 691 depending on the seed, and Bassett's own Fig 7b
  case landed on a different root in 5 of 10 identical processes. The queue is
  now seeded from the dict, which is insertion-ordered. Networks whose solve is
  basin-sensitive may now select a different (deterministic) root than they
  happened to select before; the junction scorecard is unchanged at 1700/2073
  under every seed tested, against 1700 or 1711 before. The Python suite also
  reported no intermittent `xpass` across three consecutive runs, where it
  previously did.

### Documentation
- **Corrected the `BorderCarnotLossElement` / momentum-CV accuracy claim, which
  was stated as fact in five places but is not established.** The Python
  docstring, `include/solver_interface.h`, `include/multi_port_chamber.h`,
  `docs/API_CPP.md` and the v3 addendum all asserted that
  `L = 4*(1 - cos((3/4)*delta))^2` "reproduces Hager `xi_l` and Bassett `K_inc`
  exactly at M -> 0". The repo's own Tier-1 tests
  (`python/tests/test_momentum_cv_tier1_bassett.py`) are `xfail` and contradict
  it. All five now state the anchoring as design intent with the measured
  deviation and a pointer to issue #272, so nobody takes the element as
  quantitatively validated against Hager/Bassett. **No behaviour change** --
  the physics and the coefficient are untouched; only the accuracy claim is
  corrected. Also refreshed the Tier-1 module docstring, which described a
  model version predating the sin^2(theta) + cross-coupling integration
  (including a sign-flipped `K_str`) and framed the disagreement as one
  undisentangled question; it now carries the measured K table and separates
  the two independent problems, and `experimental_bounded_solver.py` records
  that bounds do *not* explain the gap.
- **Documented the ejector's operating-regime limitations and the planned
  extension.** A new design note (`validation/ejector/OPERATING_REGIMES_DESIGN.md`)
  diagnoses a real non-converging low-primary-flow case: the critical-only
  `EjectorElement` hard-codes a choked primary in `R0` and converges to a
  self-contradictory root (primary Pt *below* the back pressure) whenever the
  forced primary flow is below the choke threshold, instead of the physical
  unchoked/subsonic branch (primary Pt just *above* the back pressure,
  `dp >= 0`). It proposes the fix -- a choked/unchoked `R0` reusing the
  compressible-orifice nozzle `(f, J)`, a C1 smooth-min critical/subcritical
  `R1` drooping omega between `P_c*` and the dead-head pressure `P_b0`, and a
  subsonic jet-pump mode -- all in one Newton-solvable residual system. No
  behaviour change; `ejector_element.py`, `_ejector_huang1999.py`,
  `include/ejector.h`, and `validation/ejector/README.md` now cross-reference
  the note and record the current limitation.
- **Identified and committed the off-design validation sources for the ejector
  regime extension.** The design note's validation strategy
  (`OPERATING_REGIMES_DESIGN.md` sec 4) now records that the subcritical-droop
  branch has air-fluid experiment + CFD data in repo (Henry et al. 2007,
  HEFAT2007; Akbarnejad & Ziabasharhagh 2025 for numeric `P_c*`/`P_b0`
  anchors), while the unchoked jet-pump branch has no dataset and is validated
  analytically (isentropic-nozzle exactness + conservation/continuity). Henry
  et al. Fig. 5 (measured air-ejector characteristic omega(Pb/P2) at 4/5/6 bar
  primary) is digitized and committed as
  `validation/ejector/data/henry2007_fig5.py`, cross-checked against the source
  figure. Corrects the earlier "Bartosiewicz 2004" citation (the source PDF is
  HEFAT2007). No behaviour change.

## [0.4.2] - 2026-07-09

### Fixed
- **`MultiPortChamberElement` degenerate iterates now carry a Jacobian, dead ports are
  classified consistently, and junctions with more than three ports are
  refused (issue #271).** Measured on the scorecard before changing anything:
  the pre-check fallback for an all-same-sign iterate returned a lossless
  residual with an *empty* Jacobian on 145 iterates from 26 solves, and none
  of those 26 converged; a port at exactly zero flow (typically an initial
  guess) was excluded from the element's supplier/collector masks but
  classified as a supplier by the Mynard closure, so its `K` came back
  mis-shaped and was silently zeroed -- a lossless junction at initialization,
  384 times per scorecard. Now: an all-ports-zero iterate returns the
  continuity residual *with* its Jacobian; an all-same-sign iterate is
  treated as the wrong-direction state it is and reaches the strict raise or
  the soft barrier, both of which pull the offending port back with a real
  derivative; an excluded port is snapped to a tiny flow in its declared
  direction before the closure sees it, so both classifications agree by
  construction and the dead-branch limit `K -> 1` applies. `MultiPortChamberElement`
  (and `ConstantKTeeElement`) raise `ValueError` at construction for more than
  three ports, since Mynard's `K` conversion is defined for three branches
  only and a 4-port junction returned finite residuals with an all-zero loss
  Jacobian. **Converged results are unchanged** -- the scorecard reproduces
  the previous tables to the digit and the 26 doomed solves remain unconverged
  (the empty Jacobian was wrong but was not what killed them).
- **`MultiPortChamberElement` no longer turns an exception inside the Mynard closure
  into a silent lossless junction (issue #271).** Three call sites wrapped
  `junction_loss_coefficient` in `except Exception`; the residual site
  returned a lossless continuity residual with an *empty* Jacobian for any
  error at all, and the FD-fallback site left a Jacobian column at zero.
  Measured across 2073 scorecard records (three sources, three topologies)
  and the full test suite, those handlers never caught a legitimate degenerate
  flow split -- the guard ahead of the call already routes those to the
  continuity fallback -- and the only thing they ever caught was a plumbing
  error, which they disguised as a plausible-looking lossless junction
  mid-solve. Now: the closure's two possible exceptions (`IndexError`,
  `ValueError`) are re-raised as a `RuntimeError` naming the element, the
  port mass flows and the supplier/collector masks, with the cause chained;
  every other exception propagates untouched. `diagnostics()` narrows the
  same way and, post-solve, reports a `closure_error` key instead of silently
  omitting the `K` fields. **Converged results are unchanged** -- the
  scorecard reproduces the previous tables to the digit -- because the
  removed path was never taken on a legitimate state.
- **`MultiPortChamberElement`'s Jacobian was missing the common-port static-pressure
  column (issue #271).** The Mynard loss term `K_i * q_dyn_com` depends on the
  common port's density, and hence on its static pressure, through two paths:
  `q_dyn_com ~ 1/rho`, and Mynard's `K` itself being a function of the port
  velocities `U = -mdot/(rho*A)`. Neither was differentiated, so the row
  carried a silently zero entry for a real solver unknown. This is the same
  defect PR #230 fixed in `ConstantKTeeElement`; it survived here because
  `test_mpce_v2_jacobian.py` FD-checks `dKQ_dmdot_separating_T` in isolation
  and nothing pinned the assembled row against the residual it differentiates.
  Taking only the first path is wrong by several per cent -- the second more
  than cancels the naive `-K*q/P` -- so the analytic form is
  `dR_i/dP_com = sign * [K_i*q_dyn/P - (mdot_com/P) * dKQ_i/dmdot_com]`, which
  reuses the existing `dKQ/dmdot` column. **Converged roots are unchanged**
  (a Jacobian correction moves the solver path, not the solution): scored
  against the Bassett separating curves the K values are identical to 3
  decimals. New whole-row FD guardrail in
  `python/tests/test_mpce_v2_jacobian_rows.py` covers both `MultiPortChamberElement` and
  `ConstantKTeeElement`.
- **CI: `publish-gui.yml` verify job had two more release-blocking races,
  both hit on the v0.4.1 release attempt.** (1) The "wait for combaero"
  step polled PyPI's JSON API, which updated several minutes before the
  Simple Repository index that `uv pip install` actually resolves
  against (separate CDN caches); the install step silently resolved the
  prior, broken `combaero` version instead of failing outright. Now
  polls `/simple/combaero/` directly. (2) The smoke test's fixed `sleep
  4` before the health-check `curl` flaked repeatedly once server
  startup (heavier now with the MPCE junction models) crept past 4s;
  replaced with a poll loop (up to 10s).

## [0.4.1] - 2026-07-08

### Fixed
- **`MultiPortChamberElement`/`ConstantKTeeElement` raised `ModuleNotFoundError: No
  module named 'validation'` on a real PyPI install.** The Mynard
  Unified0D loss-coefficient implementation lived in the dev-only
  `validation/` tree (not shipped in the wheel) and was imported
  directly by production code; a fresh `pip install combaero` (or
  `combaero-gui`, which depends on it) hit the error the moment either
  class was constructed. First surfaced by the v0.4.0 release's GUI
  publish verify job, which installs the built wheel into a clean venv
  rather than running from a repo checkout. The implementation now
  lives in `combaero.network._mynard2010`; `validation/junction/models/
  mynard2010.py` re-exports it so existing validation scripts/tests are
  unaffected. Verified by installing a freshly built wheel into an
  isolated venv outside the repo tree.
- **CI: `publish-gui.yml` no longer races `publish.yml` to PyPI.** Both
  workflows trigger off the same `v*` tag, but the core `combaero` wheel
  build (3 OSes via `cibuildwheel`) is slower than the GUI build; the
  verify job's `uv pip install` could resolve before the matching
  `combaero` version existed on the index ("No solution found" --
  first hit on the v0.4.0 release). The verify job now polls PyPI for
  the matching `combaero` version (up to 10 minutes) before installing.

## [0.4.0] - 2026-07-08

### Removed
- **GUI: the deprecated "Incomp. Warmstart" init strategy is no longer
  offered in the dropdown.** The inlet-referenced incompressible warm
  start occasionally works (heavy MFB-driven nets), but its successor
  -- the outlet-referenced variant, which tracks the compressible root
  6-7x closer and remains selectable -- covers that niche, and the
  auto-retry ladder still falls back to the inlet-referenced proxy
  internally when the outlet-referenced one starves. Legacy save
  files that stored the old strategy are migrated to
  ``outletref_warmstart`` on load. The Python API still accepts
  ``init_strategy="incompressible_warmstart"`` (with its
  DeprecationWarning) for headless workflows.

### Added
- **GUI/Network: constant-K junction model tier for the tee.** The
  tee Inspector gains a "Junction Model" selector: the default
  "Momentum CV (Mynard)" full closure, or "Constant K (handbook)" --
  fixed loss coefficients ``K_straight`` / ``K_branch`` (defaults
  0.4 / 1.0, representative for a sharp 90-deg tee) referenced to the
  common-leg dynamic head (Idelchik convention), independent of the
  flow split. Backed by the new ``ConstantKTeeElement`` (promoted from
  the certified-audit prototype), supporting both branch and merge
  topologies, with analytical Jacobian rows (including the dynamic
  head's density sensitivity to the common port's P and T) and the
  inherited MPCEv2 direction guard. Diagnostics report the fixed K
  values and the mass-flow ratio. At convergence the port total
  pressure drops reproduce ``K * q_common`` exactly (anchored in
  tests).
- **Solver: energy-consistency verification for v1
  ``MultiPortChamberBase`` -- mirror roots are demoted like
  MPCEv2's wrong-direction roots.** The v1 impulse rows are even in
  the port flows, so the sign-flipped image of any root is also an
  exact root, and the image of a physical root is energetically
  impossible (a passive junction collecting flow from low-Pt ports
  into a higher-Pt port). Post-#228 the fast cold path converges onto
  such mirror roots and reported success (three-port net: the 2.0 bar
  branch "pumping" into the 2.1 bar reservoir). The new
  ``verify_solution_consistent`` demotes them via the existing
  post-solve gate: no collector port's total pressure may exceed the
  single supplier's (1e-4 relative tolerance). Scope: separating mode
  only (exactly one supplier at the converged state) -- in joining
  mode the v1 impulse rows themselves can manufacture flow work (the
  documented deficiency that motivated MPCEv2, which merge networks
  must use), so policing energy there would reject the v1 merge
  model's own legitimate solutions. Direction cannot be checked
  instead: v1 legitimately supports runtime reversal (ejector
  regimes). End-to-end effect: the plain cold solve on ambiguous v1
  nets now self-heals -- the mirror root is demoted and the
  outlet-ref auto-retry reaches the physical root. The post-solve
  demotion message is generalized to cover both verification kinds.
- **Solver: LM-phase stall detection -- doomed primaries fail fast so
  the outlet-ref auto-retry gets the budget.** The stall detector
  (previously hybr-phase only) now also watches the LM fallback: when
  LM's best |F| plateaus (< 10% improvement over 5 s) while still far
  from tolerance (> 100x), the attempt is abandoned via
  ``SolverStallError`` instead of grinding unbounded. On the doomed
  three-port junction net the primary previously burned ~187 s flat at
  |F| ~ 4e3 before the outlet-referenced warm-start retry rescued it in
  3 s; it now fails ~10 s after LM flatlines and the full solve drops
  to ~25 s. Discrimination verified on recorded traces: fires at 5.5 s
  and 8.2 s into the two doomed LM phases, never on the 16 kg/s net
  whose LM fallback is the rescue (converges in 4.8 s). Also fixes a
  #224 side effect: the 5-tee net at 0.32 kg/s converged cold pre-#224
  (161 s of hybr grinding) but the handover sent it to an LM that
  plateaus -- it burned the full 180 s timeout and survived only via
  the auto-retry; the doomed burn now ends in seconds. Diagnostics
  record ``lm_stall: true``. Stall detection (BOTH phases, including
  the earlier hybr handover) now applies to cold attempts only: seeded
  solves (explicit x0 or a warm-start proxy seed) are the last retry
  rung and keep their full grind. Certified-audit evidence: 4/32
  outlet-ref seeded cases plateau far from tolerance for > 5 s and
  still dig out to the root (the LM-phase abort turned those slow
  successes into fast failures), and the hybr handover firing
  mid-seeded-solve turned two 20 s successes into failures (cases
  13/19) because LM cannot always finish from the same seed.
- **Solver/GUI: ``init_strategy="outletref_warmstart"`` is selectable.**
  The outlet-referenced incompressible proxy (the auto-retry's seed)
  can now be chosen directly as the primary initialization -- no cold
  attempt first -- for networks known to stall the cold path. Certified
  audit: 16/16 same-root on branch topologies, 11/16 on merge (prefer
  ``default``/``analytical_pt_prop`` for merge networks; on healthy
  branch nets the cold start remains ~4x faster, which is why this is
  not the default). Exposed in the GUI init-strategy dropdown.
- **CI: Windows extension import smoke test.** The MSVC Python extension
  build job now imports ``combaero`` and ``combaero._core`` via
  ``pytest python/tests/test_import.py``; previously the extension was
  built on Windows but never loaded, so DLL/ABI load failures could
  merge undetected.
- **GUI: frontend unit tests for topology result invalidation.**
  ``useStore.test.ts`` locks in the result-stripping semantics shipped
  with the GUI integrity fixes: structural edits (add/remove node or
  edge, new connection) drop all stored node/edge results while
  position/selection changes keep them.
- **Solver: stall-triggered handover from hybr to the LM fallback.**
  On compressible/tee networks hybr previously kept 65% of the
  wall-clock budget even when its best residual had flatlined; the LM
  fallback (which is what actually cracks these networks) only started
  afterwards. A plateau detector now ends the hybr phase early: if the
  best |F| improves by less than 10% over a trailing 5 s window while
  still far from tolerance (> 100x), a ``SolverStallError`` (a
  ``SolverTimeoutError`` subclass, so every existing handler applies)
  hands the remaining budget to the LM fallback. Verified live: a
  16 kg/s 20 bar 5-tee network drops from 103.5 s to 23.0 s
  ("hybr stalled at 1.2e4, LM |F|=1.7e-4"), while steadily-converging
  networks are untouched (the detector never fires on them). The
  far-from-tolerance gate exists because slow-but-steady grinding near
  the root is not a stall: the outlet-ref seeded retry on a 5-tee at
  0.18 kg/s crawls at |F| ~ 1e-2 before converging, and cutting it
  over to LM turned that success into a failure during live
  verification. Diagnostics record ``stall_handover: true``.
  Alternating hybr/LM was evaluated first and rejected: chaining the
  best iterate between methods poisons LM with hybr's plateau basin
  and failed on all four benchmark cases (``tmp/alternating_probe.py``,
  ``tmp/stall_detector_sim.py``, ``tmp/stall_handover_verify.py``).
- **Solver: automatic outlet-referenced incompressible warm-start retry
  for failed cold solves on compressible MPCE networks.** Cold hybr on
  compressible tee networks shows non-monotone "slow pockets" (5-tee
  GUI network: 0.32 kg/s needs 510 evals/161 s while 0.4 kg/s needs
  276/88 s), which look like hard failures at GUI timeouts. When a cold
  ``default``/``analytical_pt_prop`` solve fails, ``NetworkSolver.solve``
  now re-solves the network once in the incompressible regime with
  element densities evaluated at the DOWNSTREAM static pressure and
  warm-starts the compressible retry from that proxy solution. The
  outlet-referenced proxy tracks the compressible solution 6-7x closer
  than the classic inlet-referenced incompressible solve (RMS junction
  pressures 5-30 kPa vs 32-211 kPa across 3 diagnosed networks) and its
  seed converges cases where the inlet-referenced seed fails outright.
  The worst slow-pocket case now converges in 86 s total (vs 161 s+
  cold). Budget split: primary 40%, proxy <= 30% (up to 60 s -- heavy
  high-Re networks need ~35 s even incompressible), seeded retry the
  rest; when the stiffer outlet-referenced proxy does not converge the
  retry falls back to the classic inlet-referenced proxy seed (which
  rescued a 16 kg/s 20 bar 5-tee network where outlet-ref starved).
  The effective strategy is recorded as
  ``solver_settings_used.init_strategy = "outletref_warmstart"`` (or
  ``"inletref_warmstart"`` for the fallback) with ``auto_retry: true``.
  Opt out with ``solve(auto_retry=False)``. Evidence:
  ``tmp/outlet_ref_seed_experiment.py`` (2026-07-05) and the certified
  inverse-design audit: seeded outlet-ref scores 27/32 same-root on
  PB-driven fixtures (weak on merge topologies) vs 32/32 for the
  analytical default -- hence retry, not primary. NOTE: the homotopy
  ladder was also evaluated as the retry strategy and rejected -- its
  true wall cost is 64-200 s on the 5-tee network (each rung is a full
  solve; earlier per-rung timings undercounted it). A useful side
  effect of the 40% primary budget: hybr hands over to the
  LM fallback earlier, which alone cut the 16 kg/s network from 104 s
  to 29 s at timeout=90.

### Fixed
- **Solver: MPCE impulse rows were residual-scaled as mass flows --
  the root cause of the "doomed primary" stall class on junction
  networks.** ``_build_residual_scales`` classified all N+1
  ``MultiPortChamberBase`` rows as ``ref_mdot``, but the N impulse
  rows ``(P_i + rho_i u_i^2) - P_jct`` are PRESSURE-magnitude: in the
  scaled system hybr/LM actually optimize they were over-weighted by
  ``ref_p/ref_mdot`` (~1e5; measured 7.5e4 vs O(0.1) for correctly
  classified rows on the three-port fixture at x0). With the corrected
  classification (``[ref_p] * N + [ref_mdot]``) every known
  doomed-primary network converges COLD with no stall handover, no LM
  fallback, and no warm-start retry: three-port fixture 187 s doomed
  grind -> 2.6 s, 5-tee slow pockets at 0.26/0.32/0.40 kg/s (all
  previously failing cold) -> 5-7 s, 16 kg/s 20 bar net 23 s
  (stall+LM rescue) -> 12 s direct, one-tee choke net 25 s -> 2.3 s.
  Certified audit unchanged at 32/32 (39 s). The stall-handover and
  outlet-ref auto-retry machinery remains as a safety net but no
  longer fires on these networks. Found while investigating why the
  phi-closure prototype's trf solver crawled: least-squares objectives
  are brutally sensitive to row weights where root-finders merely
  degrade. NOTE: the faster cold path exposes that the sign-symmetric
  v1 MPCE (no ``flow_direction`` constraint, no direction verification
  hook) can converge onto exact mirror/reversed roots on
  pressure-driven nets with ambiguous direction -- a pre-existing
  model property; affected tests now pin the physical basin
  explicitly, and GUI networks are unaffected (MPCEv2 with direction
  verification).
- **GUI: editing a value could apply it to two nodes.** React Flow
  prevents default on node mousedown, so a focused Inspector input never
  blurred when clicking another node; React then reused the same input
  component (still holding the typed value) rewired to the newly
  selected node, and the eventual blur committed the old node's edit
  onto it. The Inspector's node/edge subtrees are now keyed by the
  selected entity's id, forcing a remount (fresh input state) on every
  selection change.
- **GUI: topology edits now invalidate stored results.** Adding or
  removing nodes/edges (including drag-drop creation, new connections,
  and paste) previously kept every node's ``result`` from the last
  solve, and the backend warm-starts from those -- so after adding a
  tee to a converged network the first solve was seeded with the old
  topology's solution and stalled ("poisoned by 5-tee data" on a 6-tee
  network). Structural changes now strip all stored node/edge results;
  position/selection changes keep them.
- **GUI: copy/paste no longer shares nested data between nodes.** The
  clipboard stored live node references and paste spread ``data``
  shallowly, so nested objects (composition, initial guesses, wall
  layers) stayed shared between the original and the pasted node;
  both are now deep-cloned.
- **Homotopy init: the rung ladder now shares the caller's ``timeout``
  as a total deadline.** Previously every rung received the full budget,
  so an ``init_strategy="homotopy"`` solve could run to ~6x the
  requested timeout -- blowing through the GUI's hard (soft + 10 s)
  request timeout and returning a 504 instead of a diagnosable result.

### Changed
- **Test suite and CI: 3-4x faster merge gate at unchanged coverage.**
  The four ``test_three_port_*`` network tests now share one
  module-scoped solve instead of re-running the identical ~3 min solve
  each (that single network is a known doomed-primary case: hybr
  stalls, the LM fallback grinds, and the outlet-referenced warm-start
  retry rescues it in seconds -- the shared solve pins that
  rescue-within-budget behavior with ``timeout=120``); pytest runs in
  parallel (``pytest-xdist``, ``-n auto --dist loadscope``) in CI and
  in the pre-commit hook; CI logs report ``--durations=25`` to keep
  slow-test creep visible. Python unit test job: ~29 min to a few
  minutes; local pre-commit test hook similarly.
- **CI: dependency lockfile changes now trigger the test suite.**
  ``uv.lock`` was missing from the change filters in both workflows,
  so dependency bumps (e.g. scipy/numpy) could merge without running
  a single Python test.
- **MPCE networks now default to ``analytical_pt_prop`` initialization.**
  ``NetworkSolver.solve(init_strategy="default")`` auto-upgrades to
  ``analytical_pt_prop`` when the network contains a
  ``MultiPortChamberBase``; passing any other strategy or an
  explicit ``x0`` opts out, and the strategy actually used is recorded
  in ``solver_settings_used``. Evidence: on the certified inverse-design
  audit (``tmp/mpce_audit_v2_generator.py`` + ``_runner.py``, 32
  feasible fixtures with stored ground-truth roots), analytical
  converged 32/32 to the certified root and was the fastest strategy
  overall (41 s total), vs 28/32 for the plain cold start and
  homotopy. The GUI ``SolverSettings`` schema now accepts
  ``analytical_pt_prop`` explicitly. NOTE: the earlier LHS-32 audit
  numbers (complementary default/analytical partition) are void -- that
  fixture family was later shown to be infeasible for any physical tee
  closure (equal-Pt sinks + constant-area legs demand static recovery
  the model set cannot represent).

### Deprecated
- **``init_strategy="incompressible_warmstart"``** now emits a
  ``DeprecationWarning``. On the certified audit it converged 10/32 at
  roughly 10x the wall time of ``analytical_pt_prop`` (409 s vs 41 s
  total); on a real GUI tee network it needed 445 evaluations vs 31
  for the cold start. The incompressible proxy seeds a low-residual but
  poorly conditioned region for the compressible restart. Will be
  removed once the GUI drops the option.

### Added
- **GUI: all supported init strategies selectable in the sidebar.**
  ``analytical_pt_prop`` added to the strategy dropdown (the backend
  accepted it since the auto-upgrade PR but the frontend never offered
  it); ``incompressible_warmstart`` is labelled "(deprecated)".
- **GUI: the solver status badge shows the init strategy that actually
  ran** (from ``solver_settings_used``), flagged "(auto)" when the
  solver upgraded the requested "default" -- previously the MPCE
  auto-upgrade was invisible from the GUI.

### Changed
- **GUI: "Reset All Initial Guesses" is now "Reset Guesses & Results"
  and truly resets.** It previously cleared only
  ``node.data.initial_guess`` -- but the graph builder warm-starts from
  ``node.data.result`` whenever ``initial_guess`` is empty, so the old
  button silently RE-ENABLED result-based seeding instead of giving a
  cold start. It now clears stored solve results as well; the next
  solve behaves like a freshly built network.
- **GUI: ``mpce_tee`` is displayed as "Tee".** Sidebar drag entry, node
  card default label, and Inspector header now say "Tee"; the Inspector
  section carries a short explainer naming the underlying model
  (momentum-based multi-port chamber junction, MPCE / Mynard closure).
  The internal node type and save-file format are unchanged.

### Fixed
- **GUI: failed solve results no longer poison subsequent solves.**
  The graph builder warm-started every solve from the previous result
  stored on each node -- including the best iterate of a FAILED solve
  (e.g. a timeout). Seeding from a stall point parks Newton exactly
  where the last attempt died, so one failure made a network
  "practically unsolvable" even when it converges fine from cold
  (observed on a single-tee mdot-to-pressure network: cold start 173
  evaluations, warm-from-failed stalls at |F| ~ 1.5e2 with the same
  residual signature on every retry). ``_guess_from_prior_result`` now
  ignores results with ``success != True``. Added a practical
  ``flow -> tee -> orifices -> pressure`` cold-start regression fixture
  (GUI-parity geometry) in both regimes.
- **GUI: the hard-timeout error message no longer recommends the
  homotopy init strategy** (it could appear while already running
  homotopy, and suggested a "shorter timeout" where a longer one is
  what helps).
- **GUI: inherited ambient pressure placeholder showed ``0.00``.** For
  humid-air compositions inheriting ambient conditions from the global
  solver settings, the Inspector's Ambient Pressure field displayed a
  hardcoded ``0.00`` grey placeholder while Temperature and Relative
  Humidity correctly showed the effective inherited values.
  ``UnitInput`` now accepts an effective-value placeholder (in Pa,
  scaled to the selected display unit).
- **The MPCE ``analytical_pt_prop`` auto-upgrade now requires >= 2
  ``PressureBoundary`` nodes.** The strategy's Pt propagation needs a
  boundary pressure gradient to be meaningful. On MFB-driven networks
  with a single PB (the common GUI pattern: mass inlet -> tee tree ->
  one pressure sink), propagated per-node Pt differences are
  path-length artifacts and the Bernoulli sign logic seeded REVERSED
  channel m_dots -- parking the first iterate inside the junctions'
  soft-barrier penalty region (|F(x0)| ~ 1.15e6 of pure penalty on a
  real 5-tee GUI network, hybr immediately "not making good progress")
  while the plain cold start converges in 16 evaluations. Such
  networks now keep the plain default initialization. The certified
  audit result behind the auto-upgrade is unaffected (all 32 fixtures
  are PB-driven).
- **GUI ``mpce_tee`` junctions are now built with ``strict=False``.**
  Strict mode raises whenever a Newton ITERATE (not the converged
  answer) explores a wrong-sign port m_dot, killing GUI solves that
  would otherwise converge with "Unexpected error during residual
  evaluation: MultiPortChamberElement ... declared flow_direction=... but
  observed wrong flow direction". Soft mode lets the solver explore
  wrong-direction states transiently; the post-solve
  ``verify_solution_consistent`` guard (see entry below) rejects any
  wrong-direction artifact root at convergence, so this switch does
  not admit silently wrong answers.
- **Soft-barrier artifact roots are demoted to honest failures.**
  ``MultiPortChamberElement`` with ``strict=False`` replaces a wrong-direction
  port's physics row with Pt continuity plus a one-sided penalty
  ``alpha * mdot^2`` -- and that penalty can exactly cancel the
  continuity mismatch, so the solver converges (|F| ~ 1e-10) to a state
  with a slightly reversed feed that the strict residual rejects
  (found via constant-K warm starts on the certified audit: a merge
  feed reversed to -0.0035 kg/s while the certified all-forward root
  has it at +0.213). ``NetworkSolver.solve`` now runs the element's
  ``verify_solution_consistent`` post-solve and reports such solutions
  as non-converged with an explanatory message, instead of returning a
  silently wrong flow split. This is the prerequisite for building GUI
  junctions with ``strict=False`` (Newton may explore wrong-direction
  states transiently without dying, and artifact roots can no longer
  masquerade as answers).

### Removed
- **``tee_junction`` (Bassett K-closure) retired from the GUI** (#177).
  ``MPCE Tee`` (``mpce_tee``, momentum-CV closure) is now the sole
  user-drawable tee type. Removed the Sidebar drag-create entry,
  ``TeeJunctionNode.tsx`` node card, ``TeeJunctionData`` schema,
  ``graph_builder`` dispatch, ``runner`` state projection, and all
  Inspector conditionals + display fields. The ``TeeJunctionElement``
  Python class + Python tests + validation adapters
  (``validation/junction/models/tee_junction_element_network.py``,
  ``tee_junction_raw.py``) are retained as the M -> 0 regression
  baseline for the junction validation suite -- the K-closure is no
  longer a competing user-facing model but remains authoritative for
  regression against Bassett 2001. GUI networks saved with a
  ``tee_junction`` node no longer load; they need to be recreated with
  ``mpce_tee``.

### Added
- **``init_strategy="analytical_pt_prop"`` on ``NetworkSolver.solve``**.
  Opt-in initialization strategy targeting MPCE compressible cold-start
  cases where the solver's default per-element uniform-dP propagation
  lands Newton in the wrong basin. Seeds two element unknowns the
  default path leaves at arbitrary constants:
  ``ChannelElement.m_dot`` from a Bernoulli estimate on the propagated
  Pt drop (default fallback is ``ref["m_dot"] = 0.1``), and
  ``MultiPortChamberBase.P_jct`` from the propagated Pt at the
  junction's "common" port (default fallback is ``min(boundary Pt)``,
  almost never near the true junction static). User-provided
  ``initial_guess`` entries take precedence; single-shot semantics
  (overrides cleared after ``_build_x0``). LHS-32 audit
  (``tmp/mpce_cold_start_audit.py``) partitions cases where the new
  strategy is complementary to ``default``: high-theta (> 50 deg)
  topologies favour analytical, low-theta favour default. Combined
  success rate 21/32 = 66% vs 14/32 = 44% for either alone; ship as an
  alternative attractor, not a default replacement. Follow-up work:
  fallback/auto strategy that tries multiple starts.
- **``__convergence_history__`` field on NetworkSolver solution dict**.
  Previously the per-iteration residual-norm trace lived only on
  ``solver._diagnostic_data`` (consumed by the GUI API path) and was
  invisible to plain ``solver.solve()`` callers. Now also emitted as
  ``sol["__convergence_history__"]`` so validation scripts / Python
  notebooks can plot the |F| vs iteration curve without reaching into
  solver internals. Format: list of ``{eval, t, norm}`` dicts, throttled
  at 0.1s intervals (or whenever the iterate improves), capped at 1000
  points. Both sources read from the same ``_history`` accumulator;
  this just exposes it via the public sol dict too.
- **CSV download button on ``ConvergenceChart`` GUI modal**. The chart
  already showed live the residual-norm trace + worst-residuals table
  after each solve; users can now click a download icon to save the
  raw trace as ``convergence-{timestamp}.csv`` (columns ``eval,t_s,norm``,
  plus the worst-residuals list as a comment block at the end). Lets
  you grab a slow-convergence trace, close the modal, and analyse
  later without re-solving.

### Changed
- **MPCE-v2 default Jacobian method flipped from ``"fd"`` to ``"sympy"``**
  (``MultiPortChamberElement.jacobian_method``). The sympy analytical Jacobian for
  the canonical 3-port separating T has been wired since the prototype
  PR but was kept disabled because the iteration-4 commit message
  noted FD's tiny truncation acting as regularization on some hard
  low-q mfb_two_pb cases. After the soft-barrier iteration and PR #197
  joining-side calibration shipped, that empirical concern no longer
  shows up: a focused audit comparing both methods across the full
  separating measured grid (Bassett K6 + K5, 105 points x 3 topologies)
  shows IDENTICAL convergence counts (94/32/11 of 105 across imposed_q
  / three_pb / mfb_two_pb) and IDENTICAL mean/max K errors, with
  ~1.14x wall-time speedup from cutting 4 Mynard evaluations to 1 per
  Newton step on the canonical case. Joining flow and non-canonical
  geometries (theta != 0 on straight, N != 3, port 0 not supplier)
  auto-fall-back to FD via the existing canonical-check, so the flip
  only affects separating canonical cases. Set ``MultiPortChamberElement.jacobian_method = "fd"``
  on a per-instance basis to revert.

### Added
- **K diagnostics for MPCE-v2 tee element**. ``MultiPortChamberElement.diagnostics``
  re-evaluates Mynard at the converged state and emits per-port K values
  plus topology-aware named aliases:
    - separating (``flow_direction='branch'``): ``K_straight``, ``K_branch``,
      ``mass_flow_ratio = m_dot_branch / m_dot_com``
    - joining (``flow_direction='merge'``): ``K11``, ``K12``,
      ``mass_flow_ratio`` (same convention as Bassett joining-T q)
  Plus convenience ``Pt_jct`` field (renames the legacy ``P_jct`` slot
  which already stores Pt in MPCE-v2). GUI Inspector / save files can now
  show K coefficients in the same shape that ``TeeJunctionElement``
  emits, enabling side-by-side comparison and validation against Bassett
  / Hager / Idelchik reference data without back-calculating from port
  pressures. Falls back to the parent fields only when ``port_mdots`` is
  not provided (legacy callers).
- Per-port ``port_{i}_m_dot`` field on ``MultiPortChamberBase.diagnostics``
  when the solver supplies ``port_mdots``.

### Fixed
- **``L_choke`` is interpolated within the breaking march step instead of
  snapping to the RK4 grid.** ``fanno_channel`` / ``fanno_channel_rough``
  reported the length-to-choking quantized to ``L / n_steps`` (the
  position of the last completed step). The in-channel choke barrier
  (see the barrier entry below) builds its ramp on ``L_choke``, so the
  channel's dP(m_dot) closure was a staircase near choke onset: jumps of
  ``kappa * P_in / n_steps`` (~1 kPa at 1 bar with the default 100-step
  march) every ~0.008 kg/s. Any pressure-driven network whose residual
  zero fell inside a staircase jump had NO root -- a single air channel
  (L=1 m, D=0.05 m) between two pressure boundaries failed to converge
  at EVERY Pt ratio from 1.2 to 10 under every method/init combination,
  stalling at a tread edge (best |F| ~ 90 at ratio 1.2). The sonic
  station is now interpolated: linearly in Mach between the last two
  steps for the ``M >= kFannoChokeMach`` exit, and linearly in the
  pressure probe for the ``P <= 0`` RK4-stage exits. The same
  single-channel sweep now converges at all 13 sampled ratios (1.2-10),
  landing on the monotone barrier ramp at the choke boundary as the
  barrier design intended. The choke threshold Mach (0.999) is hoisted
  to ``constexpr kFannoChokeMach`` in ``compressible.h``. NOTE: for BCs
  beyond the subsonic envelope the converged ``m_dot`` sits on the
  barrier ramp slightly above the physical choked flow -- check the
  channel's ``choked`` flag when interpreting such solutions. This also
  invalidates (again) the pre-fix LHS-32 audit numbers: the 2026-07-03
  re-run showed 25/32 all-strategy failures, most of which trace to this
  staircase because the audit's channel geometry leaves the subsonic
  envelope at Pt ratio ~1.19.
- **MultiPortChamberBase ports now report their real throughflow to the
  solver's state propagation.** ``flow_at_node`` returned 0 "for safety",
  so every port-MCN fed by the junction (all collector/outflow ports)
  received zero-mass upstream streams. Three consequences, all fixed:
  (1) the collector-port MCN closure ``Pt = P + 0.5 rho v^2`` silently
  degenerated to ``Pt = P`` -- the outflow dynamic head (~29 kPa for
  0.4 kg/s through a 5 cm port at 1 bar) vanished from the bookkeeping,
  and MPCE-v2's Mynard residual received a static pressure where its
  design expects the port stagnation pressure; (2) merging streams with
  unequal temperatures mixed with all-zero weights, so the outlet
  temperature snapped to the FIRST inlet stream's value (a 0.3 kg/s @
  400 K + 0.1 kg/s @ 300 K merge reported 400 K instead of the
  enthalpy-weighted ~375 K) -- previously masked because all junction
  tests are isothermal; (3) the degenerate closure made non-physical
  network states exactly satisfiable (dead-branch / reversed-supply roots
  at |F| ~ 1e-7). Ports now resolve their outer connecting elements'
  ``m_dot`` unknowns via an index map stashed by ``NetworkSolver`` at
  unknown-vector build time; supplier ports and single-supplier
  collectors carry their own outer flow, multi-supplier collectors the
  sum of the supplier feeds, with matching ``flow_jac_at_node`` entries
  so the closure Jacobian stays consistent. Additionally,
  ``MultiPortChamberBase.resolve_topology`` pushes inferred port areas
  onto auto-sized port-MCNs (collector ports have no upstream channel to
  inherit ``Dh`` from and previously kept the 0.1 m^2 fallback, hiding
  the face velocity even when flows were known). Converged results for
  MPCE networks shift by the collector dynamic head; junction loss
  calibrations (e.g. MPCE-v2 ``joining_etransfer_alpha``) were fitted on
  the degenerate plumbing and should be revisited against the validation
  suite. Known consequence: with honest collector stagnation pressures,
  the v1 impulse rows (``sin^2(theta)`` per-port form, Bassett
  correspondence derived for separating flow) are energetically
  inconsistent for JOINING flow -- a v1 merge junction can manufacture
  flow work and may lose its all-forward root entirely.
  ``test_step_4_adding_bypass`` migrated its junctions to
  ``MultiPortChamberElement`` (the production model, joining-capable by design);
  the v1 vs K-closure comparison test now seeds the forward basin
  explicitly since v1's sign-symmetric residuals admit the mirrored
  all-reverse root from a cold start.
- **Choked compressible channels no longer report zero resistance to the
  network solver.** ``channel_compressible_mdot_and_jacobian`` passed
  through ``fanno_channel_rough``'s truncated march: a supersonic inlet
  (M >= 1) returned ``dP = 0`` with zero ``m_dot`` sensitivity, and
  in-channel choking (``L_choke < L``) collapsed the reported drop toward
  zero as the requested velocity approached sonic. The infeasible region
  therefore looked frictionless to Newton and acted as a flat attractor
  with exact spurious network roots -- e.g. the canonical 3-port MPCE
  branch tee (pt_ratio 2.0, theta 45 deg, psi 1.0) "converged" at
  |F| ~ 3e-7 to a dead branch arm with the straight arm at M ~ 1.1 and
  zero total-pressure drop. The kernel now applies a monotone
  infeasibility barrier (``kChannelChokeBarrierKappa``,
  ``kChannelChokeBarrierBeta`` in ``solver_interface.h``): in-channel
  choking adds ``kappa * P_in * (1 - L_choke/L)`` on top of the truncated
  drop, a supersonic inlet returns ``kappa * P_in`` plus quadratic growth
  in velocity, keeping ``dP(u)`` strictly increasing across feasible,
  choked, and supersonic regions so Newton is always pushed back toward
  the feasible branch. Feasible (un-choked) results are bit-identical to
  before. The previously-spurious dead-branch state now evaluates at
  |F| > 1e4 (regression-tested); such cases honestly report
  non-convergence instead of silently returning a wrong flow split.
  Notable interplay with the port-throughflow fix (#212): the |F|~3111
  cold-start stall on audit case #4 was the flat region itself (with the
  barrier alone the default cold start converges the case in seconds),
  but the honest collector closure changes the v1 impulse merge's root
  structure so the default cold start stalls again -- while
  ``analytical_pt_prop`` now converges the case directly in ~1-2 s
  instead of the earlier ~110-eval marathon. The #210 regression test
  keeps its rescue semantics with the updated timings. The LHS-32
  audit's strategy comparison predates both fixes and should be re-run
  before designing the auto/fallback init strategy.
- **Merging/splitting MomentumChamberNode now refused at build time** (#174).
  MCN's scalar ``Pt = P + 0.5 rho v^2`` closure models a single bulk
  velocity and cannot represent merging or splitting streams within the
  chamber itself; the silent-corruption case at the time of #174 was
  ``test_step_4_adding_bypass``, where the face-convention proxy
  (channel -> Pt, orifice -> static) produced a non-physical 17 kPa
  asymmetry at a high-q merging MCN. ``FlowNetwork.validate`` now raises
  ``ValueError`` for any MCN with > 1 incoming or > 1 outgoing edges,
  pointing the user at ``MultiPortChamberBase`` (momentum-CV junction
  with one MCN per port) as the correct migration target. The
  formerly-xfailed ``test_step_4_adding_bypass`` is restored as a
  passing test on the MPCE topology and now exhibits the expected
  diameter-driven split (smaller-bore bypass < main branch). Pass-through
  MCNs (1-in / 1-out, including the GUI auto-insert pattern and MPCE
  port faces) remain valid.
- **MPCE Tee GUI shows ``m: 0.000 kg/s`` on the node card and in Live
  Telemetry** even when the junction is actually flowing. Root cause:
  ``MultiPortChamberElement.diagnostics`` only emitted per-port ``port_i_m_dot``
  values, so ``runner._build_result_objects`` (which derives
  ``ElementResult.m_dot`` from ``diag["m_dot_com"]``) fell through to
  ``0.0``. Now emits the topology-aware aliases ``m_dot_com``,
  ``m_dot_branch``, ``m_dot_straight`` (matching ``TeeJunctionElement``'s
  convention) for the canonical 3-port case. The Live Telemetry panel
  also annotates ``port_X_*`` keys with the arm letter (C/S/B) so the
  port index displayed matches the C/S/B labels on the node card.
- **Tee area inheritance** (both ``tee_junction`` and ``mpce_tee``): the
  schema docstring promised ``F_C / F_branch = None`` means "inherit from
  connected channel", but the graph_builder dispatch never implemented the
  walk. New ``_resolve_tee_areas`` helper inherits ``F_C`` from the
  common-arm channel and ``F_branch`` from the branch-arm channel,
  computes ``psi = F_C / F_branch``. Explicit user values still win;
  fallback to ``F_C = 0.01`` and stored ``psi`` when no channels are
  connected. Frontend ``NetworkCanvas`` drag-create defaults no longer
  bake in ``F_C=0.01`` / ``psi=1.0`` so the schema-level ``None`` default
  propagates and inheritance kicks in for new tees. Backwards-compatible:
  existing saved networks with explicit ``F_C: 0.01`` retain that value;
  users can clear the Inspector field to opt into inheritance.

### Added
- **Idelchik 1966 joining-T dataset** (`validation/junction/data/idelchik1966/`):
  396 raw tabulated K-coefficient cells from Section VII diagrams 7-1 through
  7-7 (converging wye with `Fs+Fb > Fc, Fs = Fc`, theta=30/45/90 deg, F_b/F_c
  in {0.1..1.0}, q in {0..1}). Two-pass independent transcription with
  formula cross-check; two confirmed Idelchik typos tagged in metadata.yaml.
  Mapping: Idelchik zeta_c.b -> K12, zeta_c.s -> K11; psi = 1/(F_b/F_c)
  matches Bassett convention. Reformatted into 36 per-(psi, theta) CSVs
  under the standard dataset naming. Used as second analytical cross-anchor
  for the MPCE-v2 joining-side calibration.
- **MPCE-v2 joining-side etransfer correction** (combaero extension to
  faithful Mynard 2010): `joining_etransfer_alpha: float = 0.2` on
  `MultiPortChamberElement`. Mynard's original etransfer collapses to zero for
  converging flow (the `(1 - flow_ratio)` factor); the asymmetric supplier
  areas in joining-T networks then leave Mynard underpredicting K_avg at
  psi > 1. The correction adds `alpha * (A_max - A_min) / (A_max + A_min)`
  to etransfer, vanishing at psi=1 (preserves the equal-area baseline) and
  growing with area asymmetry. Single parameter, calibrated against Bassett
  K11/K12 analytical + Idelchik 1966 tabulated values at psi in [1.25,
  3.33] and theta in {30, 45, 90} (analytical-only anchors, measured points
  held out as independent validation). Independent validation: imposed_q
  joining-flow audit gives mean error -10%, max error -14% vs faithful
  Mynard; three_pb and mfb_two_pb topologies essentially tied. Pass
  `joining_etransfer_alpha=0.0` for faithful Mynard, or a custom float for
  re-calibrated workflows. The default value is locked by a regression
  test (`test_mpce_v2_joining_etransfer.py`).
- **Calibration script** `scripts/calibrate_mpce_v2_etransfer.py`:
  reproduces the alpha=0.2 fit and reports stratified residuals by
  (theta, psi) for any future re-calibration.

### Changed
- `validation/junction/schema.py`: extended `Kind` Literal with
  `"tabulated"` (handbook semi-empirical values, distinct from `"calc"`
  analytical formulas and `"measured"` experimental points). Extended
  `QConvention` Literal with `"idelchik"` (equivalent to `"bassett"` for
  joining type-6).
- `validation/junction/models/mynard2010.py`: `junction_loss_coefficient`
  gains optional `joining_etransfer_alpha: float = 0.0` parameter. Default
  preserves faithful-port Mynard behavior; non-zero values enable the
  combaero joining-side correction described above.
- **Python bindings for all 12 Bassett K coefficients**: `cb._core.tee_K1`
  through `tee_K12`. Previously only K5, K6, K11, K12 were exposed; the
  rest required separate C++ tests. The new bindings let the
  `validation/junction` scorecard evaluate all 12 K against measured
  data in one pass (K1, K3, K4, K7, K8, K9, K10 now contribute real
  numbers instead of being silently unsupported). `units_data.h` and
  `docs/API_CPP.md` / `docs/API_PYTHON.md` updated to match.
- **`MultiPortChamberBase` + `BorderCarnotLossElement`** (Phase 1 of the
  momentum-CV junction). Sanctioned successor to `TeeJunctionElement` for
  n-port manifolds, ejector / merge-split-reverse regimes, and high-Mach
  operation outside the K-closure's structural `K_run ~ 0.30` ceiling
  (see `docs/junction/junction_model_v3_addendum.md` Finding 6 and the
  implementation guide `docs/junction/momentum cv implementation guide.pdf`).
  The junction owns one scalar `P_jct` and emits `N` per-port impulse-function
  residuals plus a global sum-mass residual; loss content moves to companion
  per-port `BorderCarnotLossElement` instances applying the Hager (3/4)
  effective-angle correction (reproduces Hager `xi_l` / Bassett `K_inc`
  exactly at M -> 0 on a sharp-edged lateral). Both elements coexist with
  `TeeJunctionElement` as first-class options; the user selects per junction.
  Solver change: at each port-MCN the mass-conservation row is skipped (the
  junction's sum-mass residual covers conservation; otherwise the channel and
  junction declare equal-and-opposite flows at the port, producing a zero
  row). Tier-4 cross-check against the K-closure and GUI selection dropdown
  follow in subsequent commits.
- **Solver convergence chart**: after every solve a "View convergence" link appears in
  the solver status badge (both success and failure). Clicking opens a modal with a
  log-scale `||F||` vs. time chart (no new dependency — pure SVG). The hybr phase is
  shown in orange; the LM fallback phase (if triggered) in blue with a dashed boundary
  line. A "Worst residuals" table below the chart shows the top-5 unknowns by absolute
  residual value, helping identify which node or element is hardest to converge.
- **Convergence history in API response**: `SolveResponse` now includes
  `convergence_history` (list of `{eval, t, norm}` sampled at every `||F||` improvement
  and at most every 0.1 s, capped at 1 000 points), `worst_residuals` (top 10 per-unknown
  residuals at the final iterate), `solver_settings_used`, and `lm_started_at_eval`.
  Useful for LLM-assisted diagnostics: the response is self-contained with full solver
  context.

### Changed
- **`TeeJunctionElement`**: replaced Bassett 2001 empirical K tables with the
  Unified0D compressible model (Mynard & Valen-Sendstad 2015, extended to O(M²)
  by closed-form Borda-Carnot correction). Residuals now use stagnation pressure
  (p₀) and mass-flow continuity, consistent with compressible duct elements.
  Branching tees use two K-equation residuals (loss from supplier to each collector).
  Merging tees use R_sp (equal static pressure between suppliers) + R_K,com (loss
  from mass-weighted pseudodatum to the collector), matching the LaTeX derivation
  in `docs/junction/junction_model.tex` exactly. Tee port nodes in the GUI backend
  are now `MomentumChamberNode` so static pressure P can differ from Pt at junction
  ports, allowing R_sp to be satisfied at non-zero velocity.
  Branching split closure: both collector legs now use one consistent blended
  turning-loss closure `Pt_com - Pt_leg - K_turn(theta)*q_leg - beta*x_leg*K_bc*q_ref`
  (straight leg: `theta = 0` so the turning term vanishes, leaving the Borda-Carnot
  extraction term). This removes the prior over-charge of a near-closed branch
  (`K -> 1`) and lets a low-Mach `MomentumChamberNode`-bounded branching cascade reach
  a physical all-forward, mass-conserving root. Known envelope limit: the straight run
  is a low-loss through-leg whose effective loss coefficient is capped near
  `beta*4/27 ~ 0.30` of the common dynamic head, so a large imposed straight-run
  stagnation drop (and ejector/reverse behaviour generally) is outside this closure -
  a momentum-CV junction is the sanctioned successor (see
  `docs/junction/junction_model_v3_addendum.md` Finding 6 and
  `docs/junction/momentum cv implementation guide.pdf`). The merging tee is unaffected.

### Fixed
- **Junction validation scorecard**: cells where the candidate model could
  not evaluate any record (e.g. K not yet bound, paper not supported) now
  display `-` instead of `0.000` and a misleading negative `Delta_vs_ceiling`.
  Fully-unsupported cells are also excluded from the headline aggregates.
- **`tee_junction.h` K8 transcription**: K8 (joining type 4, lateral path) had
  `c / psi` where Bassett Table 2 + Eq 34 angle correction call for `psi * c`.
  Invisible at `psi = 1` (the only point covered by the existing
  `K7-K8` complementary-identity test, since both forms reduce to the same
  expression there); diverges at `psi != 1`. Verdict surfaced by the new
  junction-validation runner (`#183`) on Bassett Fig 12b digitized data
  (`psi in {1, 2, 4}` at theta=45 deg). Corrected K8 matches Bassett's own
  published calc line in Fig 12b; the residual gap vs measured (~0.6 at
  psi=4, q=0.5) is the paper's own acknowledged moderate fit, not a code
  defect. Same fix applied to `dK8_dq`. K7-K8 identity test still passes
  unchanged. New regression test pins K8 across psi=1, 2, 4 at theta=45
  deg.
- **Channel-to-MomentumChamberNode face convention**: `ChannelElement.residuals` coupled
  the upstream node's stagnation pressure to the downstream node's *static* pressure
  (`Pt_up - P_down - dP_friction = 0`). Across an inline `MomentumChamberNode` (where
  `Pt = P + q`) this leaked the node dynamic head `q_N` as a free pressure gain - negligible
  at low Mach, fatal near choke. Since the C++ `dP_calc` is a friction-only stagnation-loss
  (the downstream static pressure is not an input), a channel now couples stagnation to
  stagnation (`Pt_up - Pt_down - dP_friction = 0`). For `PlenumNode`/boundary terminations
  (`Pt = P`) this is identical to the previous form, so only inline-MCN networks change.
  This does not address high-dynamic-head *merging* chambers, which need an axial momentum
  balance over the chamber (angle-dependent transverse loss); that closure is tracked in
  #174 and `test_step_4_adding_bypass` is `xfail` until it lands.
- **Tee residual scaling**: `_build_residual_scales` classified `TeeJunctionElement`
  rows with the mass-flow scale (`ref_mdot`, ~0.1 kg/s) although they are stagnation-
  pressure residuals (~10^5 Pa). The ~10^6 mismatch inflated the tee Jacobian rows and
  pushed the scaled condition number to ~5e8, stalling Newton. Tee rows now use the
  pressure scale (`ref_p`), dropping the scaled condition number to ~7e2.
- **Tee compressibility-correction sign**: the O(M^2) correction kappa in
  `K_dat_j_closed` (`include/tee_junction.h`) was missing the bracket factor, giving the
  wrong sign so that compressibility *raised* K. The corrected closed form
  `kappa = (s/2)[-1 - 2x/sqrt(K_inc) + s/(2 K_inc)]` is <= 0 over the physical domain
  (compressibility lowers K), matching `docs/junction/junction_model_v2.tex`.
- **LM fallback for tee networks**: the Levenberg-Marquardt fallback now also fires for
  any network containing a `TeeJunctionElement` (previously only elements flagged as
  compressible). The two supplier stagnation pressures of a tee share a near-singular
  common mode that defeats hybr; LM's regularised step resolves it, so merging tees now
  converge by default.
- **Branching-tee port areas**: a branching `TeeJunctionElement` now propagates its duct
  geometry (`F_C`, `F_branch`) to auto-area `MomentumChamberNode` ports, which cannot see
  an upstream channel during their own topology resolution and otherwise fall back to a
  sentinel area.
- **Tee phi angle convention**: the branching and merging tee models now use
  `phi_j = 0.75 * |theta_j - theta_dat|` (absolute-difference form) instead of
  the LaTeX formula `0.75*(pi - (theta_dat - theta_j))` which assumed the
  "outward-pointing" angle convention. Python passes `theta_com = theta_str = 0`
  (axis-aligned), so the old formula produced phi_str = 135° instead of 0° and
  phi_bra ≈ 168° instead of 67.5° for a 90° branch, making K values 4–8× too
  large and preventing manifold convergence. The Jacobian `dphi/d(theta_dat)` sign
  is also corrected per branch in the merging model.
- **Tee FD Jacobian eps_P scaling**: the finite-difference pressure step in
  `tee_fd_density_jacobians` is now `max(1.0, |P| × 10⁻⁴)` Pa instead of a fixed 1 Pa,
  preventing a numerically negligible step at high absolute pressures (e.g. 300 kPa inlet).
- **MomentumChamberNode Jacobian coupling for tee junctions**: `MomentumChamberNode`
  was silently dropping the `d(residual)/d(m_dot)` coupling when the upstream element
  was a `TeeJunctionElement`, because the old code emitted `{"tee.m_dot": coeff}` — a
  non-existent variable name. The solver now calls `flow_jac_at_node(nid, indices)` once
  per element in `_propagate_states()` to obtain the correct variable names
  (`tee.m_dot_com`, `tee.m_dot_branch`, or `lc.m_dot`) and propagates them via the new
  `_m_dot_jac_names` attribute on each upstream state. `MomentumChamberNode` accumulates
  these into `_upstream_m_dot_jac` and uses them in `residuals()`. This fixes a Newton
  plateau (solver stall at non-zero residual) that occurred in both incompressible and
  compressible mode whenever a `MomentumChamberNode` was connected to a tee junction.
- **Levenberg-Marquardt fallback for hybr oscillation** on compressible networks: when
  `hybr` (Powell's method) fails to converge on a compressible network — typically because
  the Jacobian is ill-conditioned near a choked operating point and Newton steps bounce
  without making progress — the solver now automatically retries from the best iterate
  using Levenberg-Marquardt.  LM's `(J^TJ + lI)` regularisation limits the step even
  when J is nearly singular, providing monotone reduction in `||F||`.  The timeout budget
  is split 65/35 between hybr and the LM fallback so both methods have working time.
  The initial hybr trust-region is also reduced (`factor=10`, down from the scipy default
  of 100) to dampen oscillating overshoots from the start.
- **Solver convergence hardening** for networks with dead-end branches (Wall nodes) and
  inherited-geometry TeeJunctions (`F_C: null`):
  - `_build_x0` now guards the TeeJunction Bernoulli estimate against `F_C=None`,
    eliminating the `TypeError: NoneType * float` crash when using *Continuation* init
    on a network whose tee junctions inherit area from adjacent channels.
  - Initial guess for dead-end channels (connected to a `WallNode`) is now seeded at
    `m_dot=0` instead of the full upstream flow, so Newton starts from the physically
    correct state.
  - `WallNode` initial pressure is now set equal to the upstream node pressure (no
    friction drop at zero flow) instead of the BFS-propagated underestimate.
  - Default `maxfev`/`maxiter` limit is now `max(500, 200*(n+1))` (matching scipy's
    hybr internal default) computed from the actual problem size, replacing the
    hard-coded 500 that caused convergence failures on larger networks or near-choking
    operating points.
  - `ChannelElement.diagnostics` applies the same Re floor (`sqrt(Re²+64²)`) used by
    the C++ friction correlations, preventing `f` from reporting astronomically large
    values (e.g. 5.7 × 10¹⁶) for dead-end channels with near-zero mass flow.

### Added
- **Copy-paste nodes** (Ctrl/Cmd+C / Ctrl/Cmd+V): selected nodes can be duplicated
  on the canvas. Physics and configuration data (including orientation) are copied;
  custom labels and stale solver results are not. Connections are not duplicated —
  the pasted nodes start unconnected. Pasted nodes land 50 px offset from the
  originals and are immediately selected so they can be repositioned.
- **Upstream geometry inheritance**: geometry fields on elements and nodes can now
  be left blank (schema value `null`); the GUI Inspector shows an inherited value
  as a greyed placeholder with a "Reset to inherited" / "Reset to circular" button.
  Affected fields:
  - `OrificeElement.diameter` — inherits from adjacent `ChannelElement`
  - `MomentumChamberNode.Dh` — inherits from upstream `ChannelElement.D`
  - `CombustorNode.Dh` — inherits from upstream `ChannelElement.D`
  - `AreaChangeElement.F0` / `F1` — inherit from upstream/downstream `ChannelElement`
  - `TeeJunctionElement.F_C` — inherits from common- or straight-arm `ChannelElement`
  - `TeeJunctionElement.F_branch` — new field; inherits from branch-arm `ChannelElement`;
    `psi = F_C / F_branch` is now computed automatically from the two areas, replacing
    the manual psi input
  - `ChannelElement.D` — inherits from upstream `AreaChangeElement.F1`,
    `ChannelElement`, or `TeeJunctionElement` (correct arm area used for branch vs. straight)
  - `ChannelElement.Dh` — new optional field; defaults to `D` (circular assumption);
    override for non-circular ducts to use correct hydraulic diameter in friction/Nu
- **TeeJunction Inspector overhaul**: three-arm area layout (C, S read-only, B) with
  per-arm inheritance and a derived ratio card showing `psi`, `D_C → D_B`, and an
  expansion/contraction label — mirrors the `AreaChange` inspector style. The `F_straight
  = F_C` Bassett (2001) model constraint is now clearly annotated.
- **AreaChange Inspector**: non-editable area ratio summary card showing `F₁/F₀`,
  expansion/contraction/straight label, and equivalent circular diameters `D₀ → D₁`.
- **Discrete loss area inheritance**: `PressureLossElement` now inherits flow
  area from an adjacent `ChannelElement` when none is explicitly set, so discrete
  losses can be placed after any element (not only combustors / momentum chambers).
  The Inspector placeholder and hint text update automatically when a channel is
  connected upstream.
- **Topology resolution robustness**: two-pass topology resolution now handles
  mutual dependencies (e.g. channel and tee both with null geometry) without
  ghost fallback values; deferred fallbacks applied only after second pass.
- `WallNode`: closed-end (zero mass-flow) boundary node for modelling dead-end
  branches and symmetry planes of ring manifolds. Drag from the GUI palette or
  instantiate directly via `combaero.network.WallNode("id")`.

## [0.3.1] - 2026-05-21

### Fixed
- `combaero-gui` dependency pin updated from `~=0.2` to `~=0.3`; the old pin
  resolved to `combaero 0.2.6` on PyPI which predates `TeeJunctionElement` and
  `VortexElement`, causing an `ImportError` at startup (#156).
- `publish-gui.yml` CI workflow upgraded from Node 20 to Node 22; `pnpm latest`
  (v11) requires Node ≥22 and crashed with `ERR_UNKNOWN_BUILTIN_MODULE: node:sqlite`
  on Node 20 (#157).

## [0.3.0] - 2026-05-20

### Changed
- **Solver hardening** against unphysical intermediate states during Newton iteration:
  - Temperature floor raised from 50 K to 200 K (both in `_propagate_states` and
    `_get_node_state`); ceiling of 5000 K added to bound combustion runaway.
  - Static pressure floor of 1000 Pa added in `_get_node_state` to prevent NaN
    propagation in thermo calls when P drifts near zero.
  - Penalty residual (triggered on unphysical state exceptions) now encodes
    `F = x_scaled − x_best_scaled` with `J = I`, so the Newton step jumps back
    exactly to the best physical iterate seen so far — previously the constant
    penalty `[10,…,10]` could drive iterates further into unphysical territory.
  - `CombustorNode.compute_derived_state` now wraps the C++ combustion call in a
    try/except and re-raises with a diagnostic message; the solver penalty path
    already handled C++ exceptions, but the new wrapper provides traceability.
  - `graph_builder` now seeds `initial_guess` from the previous solve result
    (`data.result` in node/element data) when no explicit user override exists,
    giving repeated solves a warm start without requiring the `continuation`
    strategy.

### Added
- Global ambient conditions (T, P, RH) in the sidebar "Global Settings" panel, defaulting
  to ISO 2314 gas-turbine reference conditions (288.15 K, 101325 Pa, RH 0.6). When set,
  these override the per-node `ambient_T`, `ambient_P`, and `relative_humidity` fields on
  all boundary nodes that use the `humid_air` composition source. The CompositionEditor
  shows a live `global: —` / `global: value` hint below each affected field.
- All new boundary nodes now default to `humid_air` composition (was `dry_air`). Existing
  saved networks with explicit `source: "dry_air"` are unaffected.
- `test_humid_air_rh0_equals_dry_air`: verifies `species.humid_air_mass(T, P, RH=0)` is
  identical to `species.dry_air_mass()` within floating-point tolerance.
- Global Nu/f multipliers in the sidebar "Global Settings" panel. Setting either
  scales all channel, combustor, momentum-chamber, and discrete-loss elements
  multiplicatively on top of their per-element values — useful for network-wide
  sensitivity studies without touching individual elements.
- `VortexElement` network element: models centrifugal pressure rise in rotating cavities
  (disc pumps, pre-swirl chambers, inter-stage labyrinth spaces) using the Vatistas n-vortex
  model. Parameters: `r_c` (core radius), `r_out`, `r_in` (evaluation radii), `omega_rpm`
  (shaft speed; `None` = inherit global setting), `n` (Vatistas shape parameter, default 2).
  Residual `Pt_out - Pt_in - dP_vortex = 0` with analytical Jacobians w.r.t. Pt and inlet
  density. GUI: Tornado icon (violet), inspectable from the sidebar with a global shaft-speed
  field and per-element override.
- Vatistas (1991) n-vortex model: `vatistas_v0_bar`, `vatistas_vr_bar`,
  `vatistas_pressure_integral` and their analytical derivatives, plus dimensional
  `vatistas_v_theta` / `vatistas_delta_p` with full analytical Jacobians w.r.t.
  r, Γ and r_c. Closed-form antiderivatives for n=1 and n=2; composite Simpson
  quadrature for general n. Python class `VatistasVortex` wraps the free functions
  with vectorised `V_theta`, `dV_theta_dr`, `delta_P`, and `delta_P_bar` methods.
  Validated against digitised experimental data from the original paper
  (Vatistas, Kozel, Mih, *Exp. Fluids* 11, 73–76, 1991).
- Channel element friction model selector in the GUI Inspector: Haaland (default),
  Serghides, Colebrook-White, or Petukhov. The Roughness field is greyed out when
  Petukhov is selected (smooth-pipe model, ignores roughness). Petukhov is the natural
  pairing for channels that already use the Petukhov/Gnielinski heat transfer correlation.
- `NetworkRunner` and `NetworkResult` in `gui/backend/runner.py`: programmatic
  network driver for use in Python scripts and Jupyter notebooks without
  starting the web server.
  - `NetworkRunner.from_file(path)` / `from_dict(d)`: load a GUI-saved JSON.
  - `NetworkRunner.solve(overrides, method, init_strategy, timeout)`: stateless
    single solve with optional boundary-condition overrides keyed by node label
    (``"air_inlet.m_dot": 1.2``).  Each call deep-copies the schema so repeated
    calls inside optimisation loops do not interfere.
  - `NetworkRunner.sweep(params, metrics, ...)`: parametric sweep over a
    DataFrame of BC overrides; returns a compact result DataFrame (one row per
    solve) when ``metrics`` is given, or the full per-entity detail export
    stacked with a ``_sweep_index`` column when omitted.
  - `NetworkResult.get(key)`: scalar result by raw solver key **or by
    ``<label>.<quantity>`` format** (e.g. ``result.get("combustor.T")``),
    resolving GUI node labels to internal IDs automatically.
  - `NetworkResult.node_state(label)`: full solved thermodynamic state dict for
    a node by GUI label.
  - `NetworkResult.to_dataframe()`: full-detail DataFrame identical to the GUI
    CSV export, with unit-annotated column names.
  - `NetworkResult.swap_boundary(label, new_type)`: retype a boundary node
    (mass ↔ pressure) and auto-populate its required field from the solved
    state (``Pt`` when going mass→pressure; ``m_dot`` when going
    pressure→mass).  Returns a new ``NetworkRunner`` pre-seeded with the full
    solved state — including junction-node pressures and element flows — so the
    first re-solve starts from the exact operating point and converges without
    requiring a warm-start strategy.  Enables off-design pressure sensitivity
    studies and matched boundary-condition workflows.
- `python/examples/network_runner_bc_swap.py`: example loading the compressible
  CH4/air combustor, solving the design point at fixed air mass flow, swapping
  the air inlet to a pressure boundary, and sweeping inlet total pressure
  ±10 % around the design point (11/11 points converge).
- `_schema_maps` and `_build_result_objects` helpers extracted from the HTTP
  solve/export handlers into `runner.py`; shared by `NetworkRunner` and
  `main.py` to eliminate duplication.
- `python/examples/network_runner_combustion_sweep.py`: end-to-end example
  loading a GUI-saved compressible CH4/air combustor network, injecting node
  labels, and sweeping equivalence ratio 0.3–0.9 via `NetworkRunner.sweep()`.
- `scripts/check-python-style.sh` now scans `gui/backend/` alongside the core
  Python source directories so backend style is checked locally, not only in CI.
- `scripts/check-gui-style.sh` now supports a `--fix` flag; default (no flag)
  is read-only, matching the CI `lint-gui` job.

### Fixed
- Network solver no longer crashes on combustor networks where Newton steps temporarily
  produce non-physical states (T ≤ 0 or P ≤ 0). `adiabatic_T_complete_and_jacobian_T()`
  and `adiabatic_T_equilibrium_and_jacobians()` now clamp T and P to safe floors instead of
  throwing; `_propagate_states` and `_get_node_state` in the Python solver clamp T ≥ 50 K
  so that reverse-flow mixing artefacts (negative T_mix from negative m_dot streams) never
  propagate into physics calls.
- `friction_petukhov()` and `friction_and_jacobian_petukhov()` no longer throw
  for Re < 3000; the Petukhov formula is mathematically continuous so values below
  the validated range are returned as extrapolated rather than raising an error.
  `CorrelationValidity::EXTRAPOLATED` is set when Re < 3000; `VALID` otherwise.
- `NetworkGraphSchema` now normalises the camelCase ``solverSettings`` key
  (written by the GUI's JSON export) to ``solver_settings`` via a Pydantic
  ``model_validator``, preventing solver regime, method, and timeout settings
  from being silently dropped when loading a GUI-saved file.

### Changed
- `gui/backend/main.py`: result-building and DataFrame construction replaced
  with calls to the shared helpers in `runner.py` (~135 lines removed).

- Thermal wall GUI: arrow direction, probe hot/cold labels, dropdown options, depth hint
  ("0 = Hot, L = Cold"), and "Hot Side Wall Temp" now all reflect the actual heat flow
  direction from the solver (Q sign) rather than the topology of the drawn edge. Wiring a
  thermal wall in either direction gives correct visual feedback without manual adjustment.
- Polar molecule transport model (Monchick-Mason Omega*(2,2) table, 37 T* x 8 delta* points):
  - `Transport_Props` struct gains `dipole_moment` [Debye] and `z_rot` fields.
  - `omega22()` refactored to accept `(T_star, delta_star)` and bilinearly interpolate
    the full Monchick-Mason table (same table as Cantera's `MMCollisionInt`).
  - `compute_delta_star()` helper computes reduced dipole moment in SI.
  - Thermal conductivity upgraded from simple Eucken to Mason-Monchick modified Eucken
    with Parker Z_rot temperature correction (matches Cantera's `fitProperties` formula).
  - All 15 species entries updated with `dipole_moment` and `z_rot` (GRI-Mech 3.0 values).
  - H2O dipole was previously missing; now correctly set to 1.844 D.
  - NH3 polarizability corrected from 0.0 to 2.100 A^3.
  - `generate_thermo_data.py` updated to emit both new fields; `TransportProps` dataclass
    extended with `dipole_moment` field.
- Polar species Cantera validation tests (`TestPolarSpeciesTransport`): pure NH3 and H2O
  viscosity/conductivity vs GRI-Mech at 300/1000/2000 K; NH3/air and H2O/N2 mixtures.

### Fixed
- Compressible solver convergence with mass-flow boundary inlets: `_propagate_mdot_guess`
  was applying a Mach-0.3 cap even to elements directly seeded by a `MassFlowBoundary`,
  causing the initial guess to underestimate the known mass flow and leading hybr to
  converge to a spurious local minimum. The cap is now skipped for elements whose
  immediate upstream node is a `MassFlowBoundary`.
- Edge connection paths for rotated nodes: `FlowEdge` and `ThermalEdge` now apply a
  `rotatePosition` helper that maps each handle's static `position` prop (e.g.
  `Position.Left`) through the connected node's CSS rotation in 90° steps. This gives
  `getBezierPath` the correct perpendicular exit/entry direction after rotation — for
  example, a 180°-rotated node's `flow-target` (logically Left) is correctly treated as
  Right, producing a smooth C-curve instead of an initial vertical leg routing behind the
  element bounding box.
- Rotated-node handle placement and connector triangle orientation: each `Handle` now
  receives an inline style from `handleStyle(basePos, rotation)` in `nodeUtils.ts`.
  This overrides ReactFlow's position CSS classes (which apply in the rotated coordinate
  system and would land a handle on the wrong visual edge) and counter-rotates the handle
  by `−rotation` so the `::after` CSS triangle defined in `index.css` always points
  perpendicular to the edge it sits on.  The centering uses `calc(50% − 6px)` rather than
  `transform:translate`, so `getBoundingClientRect()` returns the correct screen centre
  regardless of node rotation — keeping edge spline endpoints accurate.
- `combaero-gui` white page on fresh PyPI install: `frontend/dist/assets/` was not bundled in
  the wheel because the `**` glob in `package-data` is unreliable below setuptools 69; added an
  explicit `assets/*` pattern and a `MANIFEST.in` as belt-and-suspenders.
- Silent NH3 composition drop in `TestTransportProperties` and `TestCombustionEquilibrium`
  Cantera validation tests: hardcoded 14-item species lists were missing NH3 (index 14), so
  any NH3 mole fraction was silently discarded when building the Cantera composition dict.
- Floating-point precision failure in `test_channel_ribbed_jacobians` on Linux: the
  `dq/dT_hot == -h` identity is exact analytically but accumulated ~4 × 10⁻¹⁴ rounding
  error after NH3 extended `dry_air()` from 14 to 15 elements; replaced exact equality with
  `abs(...) < 1e-9` (matching the pattern used elsewhere in the same file).

### Added
- Ammonia (NH3) species to the thermodynamic and transport database.
- Dynamic species count support in C++ and Python test suites to prevent regressions when modifying the database.
- Cantera validation tests fully migrated to dynamic species API (`num_species()` /
  `species_name(i)` from `combaero._core`); hardcoded species lists removed from all
  `set_cantera_composition` helpers and `SPECIES_ORDER` / `CB_SPECIES` constants.
- `scripts/test-pypi-wheel.sh`: local pre-publish wheel verification — builds both wheels,
  inspects the GUI wheel for bundled assets, installs in an isolated venv, and smoke-tests that
  the server starts and serves JS with the correct MIME type.
- CI `verify` job in `publish-gui.yml` gates the PyPI publish step on the same wheel checks.
- `test_gui_export_dataframe.py`: unit tests for the CSV export pipeline's DataFrame
  construction, column-rename, and `to_csv` logic; covers pandas 3.0 compatibility.

### Changed
- Narrowed `requires-python` to `>=3.12`; dropped `from __future__ import annotations` throughout
- Self-referential `classmethod` return types migrated to `typing.Self`

### Added
- Tee junction pressure-loss model (Bassett 2001) in `include/tee_junction.h`:
  - Pure K-coefficient functions K5, K6, K11, K12 with analytical dK/dq derivatives.
  - `merging_tee_K_straight/branch` and `branching_tee_K_straight/branch`: blended,
    smooth K functions valid for any real inputs (soft-lower protection, tanh blend for
    topology reversal, no NaN or throws outside the validated range).
  - `tee_check_inputs(q, psi, theta)` returns `TeeInputStatus` with per-parameter
    validity flags and an overall `CorrelationStatus` (Valid / Extrapolated).
  - `soft_lower(x, lo)`: smooth max(x, lo) using sqrt regularisation.
- `TeeJunctionResult` struct and `merging_tee_residuals_and_jacobian` /
  `branching_tee_residuals_and_jacobian` in `solver_interface.h/.cpp`: full Jacobians
  w.r.t. mass-flow rates; FD Jacobians w.r.t. P, T, and species mass fractions.
- Python bindings in `_core.cpp` exposing `TeeJunctionResult`, both residual functions,
  and utility functions (`tee_K5`, `tee_K6`, `tee_K11`, `tee_K12`, `tee_blend_weight`,
  `tee_check_inputs`, `merging/branching_tee_K_straight/branch`).
- GTest suite `tests/test_tee_junction.cpp` (13 tests) and Python test suite
  `python/tests/test_tee_junction.py` (19 tests) covering reference values, derivative
  accuracy, zero-flow regularisation, robustness outside the validated range, and
  Jacobian finite-difference verification.
- `TeeJunctionElement` network component: three-port flow element using Bassett 2001
  coefficients, solvable via `NetworkSolver`.  Supports both merging and branching
  tee configurations.  Unknowns are `{id}.m_dot_com` (total common-arm flow) and
  `{id}.m_dot_branch`; straight flow is implicit.  Analytical Jacobian provided for
  all unknowns including pressures, temperature, and species.
- Multi-port element protocol on `NetworkElement`: `all_source_nodes()`,
  `all_sink_nodes()`, `flow_at_node()`, `flow_jac_at_node()` extension points allow
  N-port elements to integrate with the existing solver and graph topology without
  breaking 2-port elements.
- `NetworkSolver` Bernoulli-based initial guess for `TeeJunctionElement` mass flows
  (A * sqrt(2*rho*dP)) to ensure Newton converges from a physically reasonable start.
- Network-level tee tests in `python/tests/test_tee_network.py` (15 tests): graph
  connectivity, solver convergence and mass conservation for both merging and
  branching configurations, Jacobian FD verification, and `validate()` error paths.
- `TeeJunctionElement.diagnostics()` extended to expose C++-computed `K_straight`,
  `K_branch`, and `correlation_extrapolated` (1.0 when Bassett inputs are outside the
  validated range) alongside the existing `m_dot_com/straight/branch` and `q` fields.
- `test_tee_network.py`: conservation tests for mixed-composition merging streams
  (`test_merging_tee_node_mass_conservation`, `test_merging_tee_species_conservation`,
  `test_merging_tee_enthalpy_conservation`) and full central-difference Jacobian coverage
  for both merging and branching tees including mixed-composition sensitivities.
- GUI tee junction inspector, palette item, and React Flow node visual.
- GUI solver status panel now shows full error messages in a scrollable area (previously
  truncated to a single line).
- Tee junction node now shows port labels S / C / B (straight / common / branch) in
  blue (input) or amber (output) based on the current tee type, so wiring direction
  is unambiguous at a glance.
- Inspector shows a port guide card that updates with the selected tee type and
  describes which arms are inlets vs. outlets.
- `graph_builder` auto-inserts a `PlenumNode` junction when an element is connected
  directly to a tee port, removing the requirement to manually place a plenum on every
  tee arm.

### Changed
- `TeeJunctionElement.diagnostics()` key `q` renamed to `mass_flow_ratio`
  (unit kg/kg) for clarity; updated in Python, GUI inspector, CSV export, and
  `quantities.ts` catalogue.
- `TeeJunctionElement.validate()` now accepts negative branch angles: `|theta|`
  must be in `(0, pi/2]` and `abs(theta)` is passed to the C++ correlation,
  since the model uses `cos(0.75*theta)` which is symmetric. The inspector hint
  text and frontend validation both reflect the `(-90, 0) U (0, 90]` deg domain.
- GUI FLOW DISTRIBUTION section in the tee inspector now matches the LIVE
  TELEMETRY visual style (card background, `text-xs` heading, stacked label+value
  cells with `gap-y-3`).
- `TeeJunctionData` schema gains `initial_guess: dict[str, float]` field;
  `graph_builder` wires it through `_expand_initial_guess`, and the inspector
  shows an `InitialGuessEditor` for `m_dot_com` and `m_dot_branch`.
- Frontend validation (`validation.ts`) covers tee junction `theta`, `F_C`, and
  `psi` with the same yellow warning-box UX as other elements.
- CSV export mole fractions X for tee junction common-arm state (consistent with
  other element types).

### Fixed
- `NetworkSolver` Bernoulli initial guess for `TeeJunctionElement` now uses the
  propagated `p_guess` dict for non-`PressureBoundary` nodes (plenums, mass-flow
  boundaries, momentum chambers).  Previously `getattr(node, "Pt", ref["P"])` fell back
  to the reference pressure for interior nodes whose `Pt` is a solver unknown, giving
  `dP_total = 0` and a degenerate `m_dot` starting point.
- CSV export for tee junction rows now includes the common-arm thermodynamic state
  (T, P, Pt, Tt, rho, and full species composition Y[N2] ... Y[NH3]) so that the
  canonical mixed state is accessible without joining to node rows.
- `_propagate_mdot_guess` overwrote rather than accumulated flow for multi-source elements
  (merging tee): the branch-arm visit (0.1 kg/s) silently clobbered the straight-arm visit
  (1.0 kg/s), leaving the tee initial-guess at 0.1 instead of 1.1 kg/s.  The resulting 7x
  underestimate prevented Newton from converging for small `F_C` (e.g. 0.01) where the
  quadratic pressure-drop scaling amplifies initial-guess errors.

## [0.2.6] - 2026-05-03

### Fixed
- `combaero-gui` wheel did not bundle `frontend/dist/assets/`; the `**` glob in
  `package-data` is unreliable below setuptools 69 — added explicit `assets/*`
  pattern and `MANIFEST.in` as belt-and-suspenders (#124).
- `setuptools_scm` appended a local version identifier (`+g…`) to the GUI wheel
  version, causing PyPI to reject uploads; suppressed with
  `local_scheme = "no-local-version"` (#125).
- Tracked `gui/frontend/dist/` in git caused dirty working tree during CI builds,
  making the version suffix appear; removed from index (#125).
- `publish-gui.yml` skipped after re-pushing a tag because the `workflow_run`
  trigger resolved to `failure`; replaced with a direct `push: tags: v*` trigger
  (#126).

### Added
- `scripts/test-pypi-wheel.sh`: local pre-publish wheel verification — builds both
  wheels, inspects bundled assets, installs in an isolated venv, and smoke-tests
  MIME types (#124).
- CI `verify` job in `publish-gui.yml` gates the PyPI publish step on the same
  checks (#124).

## [0.2.4] - 2026-05-02

### Fixed
- GUI bundle not rebuilt before release; SPA catch-all route incorrectly served
  JSON for unmatched paths; development banner visible in production build (#118).

## [0.2.3] - 2026-05-02

### Fixed
- Relative API URLs broke on non-root deployments (#117).
- SPA catch-all route returning JSON instead of HTML for deep-linked paths (#117).
- PyPI documentation links pointed to wrong locations (#117).

## [0.2.2] - 2026-05-02

### Fixed
- `combaero-gui` frontend assets not served when installed from PyPI; broken
  documentation code examples (#116).

## [0.2.1] - 2026-05-02

### Fixed
- `gui` extra missing from `combaero-gui` package metadata (#115).
- `gui/frontend/dist/` artefacts incorrectly tracked in git (#115).

## [0.2.0] - 2026-04-26

### Added
- Network GUI: React Flow visual network builder backed by FastAPI (`combaero-gui`)
- Analytical combustor Jacobians and full network Jacobian integration
- Continuation solver with GUI method-selection sidebar
- Stagnation and transport property helpers with C++ performance optimisations
- Inverse solver API: `calc_T_from_u`, `dg_over_RT_dT` bindings
- Multi-stream axial staging combustor network example
- Diagnostic telemetry schema (Pt/Tt stagnation properties)
- Solver hard timeouts and empty-network guard
- Area-change element diagnostics and validation hardening

### Changed
- uv workspace consolidation; project-wide ruff/uv modernisation
- pybind11 requirement raised to `>=3.0.4`
- Frontend toolchain standardised on pnpm, Vitest, Node 24 LTS
- CI: separated `gui` path filter from C++ source matrix; added `lint-gui` job

### Fixed
- State property units mismatch in example runner
- Solver stability under degenerate initial conditions

[Unreleased]: https://github.com/thirox-lab/combaero/compare/v0.11.1...HEAD
[0.11.1]: https://github.com/thirox-lab/combaero/compare/v0.11.0...v0.11.1
[0.11.0]: https://github.com/thirox-lab/combaero/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/thirox-lab/combaero/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/thirox-lab/combaero/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/thirox-lab/combaero/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/thirox-lab/combaero/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/thirox-lab/combaero/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/thirox-lab/combaero/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/thirox-lab/combaero/compare/v0.4.2...v0.5.0
[0.4.2]: https://github.com/thirox-lab/combaero/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/thirox-lab/combaero/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/thirox-lab/combaero/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com/thirox-lab/combaero/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/thirox-lab/combaero/compare/v0.2.6...v0.3.0
[0.2.6]: https://github.com/thirox-lab/combaero/compare/v0.2.4...v0.2.6
[0.2.4]: https://github.com/thirox-lab/combaero/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/thirox-lab/combaero/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/thirox-lab/combaero/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/thirox-lab/combaero/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/thirox-lab/combaero/releases/tag/v0.2.0
