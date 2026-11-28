# Exact composition tangents of the Fanno kernels (#1311, #1310): provenance record

0.13.0, theme 2 ("exact composition tangents") of the release plan (#1341). The rule this serves is CLAUDE.md's: no finite difference hidden in a kernel.

## The problem

A duct's flux G solves an implicit relation, and its composition column dG/dv (v a mass-fraction direction) comes from implicit differentiation.

- **Adiabatic** (`raw_flow`): psi(Me) = Lam(G(Me), Me) with G = P / phi(Me). Its parts are:
  - the exit pressure per unit flux phi(Me, Tt0, X);
  - the isentropic inlet flux(Pt0, Tt0, M_in, X);
  - the length between the two Mach numbers at fixed limits, ell(G, Tt0, X).
- **Heated or cooled** (`rayleigh_raw_flow`): E(G; X) = P_exit after the fixed 32-step RK4 march in theta = M^2 - M^4/2.

At low drive these sums cancel to about drive / P, and twice for a composition direction: inside psi_v, then in G_v - G_M psi_v / psi_M. That is 5e-6 at 1 Pa on 2 bar.

`d_ddir` differenced each piece along v with a 1e-5 step. That carried about 1e-8 relative error, and the cancellation amplified it to:

- **Adiabatic (#1311):** 1-8 % near the 1 Pa floor. The sub-floor blend inherited it, up to 1.93x.
- **Fanno-Rayleigh (#1310):** up to 12 % along an O2/N2 trade at 1 Pa, and 21 % along pure O2 at 1.05 Pa. The 2.8e-3 in the issue was for another duct.

A larger step does not help: it trades round-off for truncation, and #1311's 100x test showed it.

## Formulation

Forward mode along v, with no step. In `src/fanno_mach.cpp`: `MixTangent`, the `*_t` templates and `rayleigh_relation_ddir`.

- **The direction.** v moves the mole fractions by dX_i = MW (v_i / W_i - X_i sum_k v_k / W_k), the step `mass_to_mole(Y + t v)` takes. It keeps sum X = 1.
- **Linear mixture properties.** Molar cp, h, dcp/dT, d2cp/dT2 and MW are linear in X, so each one's derivative along dX is the same sum with dX in place of X. Mass-based properties carry -q MW_v / MW. gamma, a^2 and its T-slope follow from those by the chain rule.
- **Entropy.** The isentropic ratios use (s(T1) - s(T2)) / R at one composition. The mixing terms cancel there, so the derivative is sum_k dX_k (s_k(T1) - s_k(T2)) / R_u, which has no ln X term and holds for an absent species.
- **Viscosity.** `viscosity_dX` (Wilke) contracted with dX.
- **The static temperature.** Implicit in h(T) + M^2 a^2(T) / 2 = h(Tt).
- **The adiabatic length.** The same fixed 4-panel Gauss-Legendre rule applied to the integrand's tangent, which is the exact tangent of the rule's sum. The length's G and Tt0 slopes, central differences until now, take the same pass with seeds (G, Tt0).
- **The Fanno-Rayleigh march.** The same RK4 grid with a composition slot, one pass per direction, for:
  - the unchoked exit relation;
  - the exit-choked relation (the backward march from sonic less the inlet theta);
  - the inlet-choked sonic flux.
- **The sub-floor low-Mach closure.** Composition enters only through R (both prefactors) and each node's mu, so its column is closed form (`low_mach_drop_ddir`).

The existing DualN<4> state tangents of the Rayleigh march (#1276) are untouched, so the Rayleigh kernel's G and state tangents are bit-identical.

## Ground truth

- **Reference.** A two-level Richardson ladder of G itself along v, with Y steps 8e-2 to 1e-2. Near the floor, G's root carries about 2e-9 relative jitter and dlnG/dv is only about 1.5e-3, so no difference of G resolves dG/dv much below 3e-5 there.
- **Adiabatic, Pt0 2e5 Pa, Tt0 400 K, L 0.5 m, D 0.01 m, haaland, O2/N2 and pure-O2 directions, static and total exits, drives 1 to 300 Pa.** The branch is within 1.35e-4 of the reference, and within its spread everywhere. Main is off by up to 8 %.
- **Fanno-Rayleigh, the same duct heated 400 -> 500 K.** The branch is within 3e-4 (reference spread 3e-4). Main is off by up to 12 % (O2/N2) and 21 % (pure O2, 1.05 Pa).
- **At 300 Pa** both kernels agree to 2e-8.
- **Pure O2 into the heated duct near 1 Pa is a near zero crossing:** 4.6e-7 against 4.5e-4 for N2. No difference resolves it there.
- **Linearity:** (e_O2 - e_N2) / 2 equals the O2/N2 trade to 1e-8.

## Choices not forced by the spec

- **One pass per direction** with a single composition slot, rather than one wide dual over all directions. The direction count is a runtime size, and one pass costs about what d_ddir's two evaluations did. Measured: 13.4 ms against 13.8 ms for 10 directions (adiabatic) and 43 ms against 48 ms (Fanno-Rayleigh), on a loaded machine.
- **The tangent of the discrete method.** For the length and the march, these are the exact tangents of the fixed rules, not of the continuous integrals: they are the derivative of the G the kernel returns.

## Not changed

- **#1317.** The sub-floor sigma slopes still difference (`compute_floor_sigma_derivs`), the owner's exception from #1303. sigma is built from dG/dP_target, so its slope is a mixed second derivative, which needs a second-order type (nested duals). DualN holds doubles, so it does not nest as it stands.
- **The throat-choked Rayleigh branch** still differences its re-solved choked flux G*. That relation goes as (G* - G)^(2/3), which is not differentiable at the choke (#1276).
- **`d_ddir`** stays for those two callers.

## Dead ends

- **#1324's polish step in `static_T`** was the first candidate for low-drive tangent noise. It does not change the cold-vs-warm spread: the spread is the G root's own conditioning, recorded on #1324.

## Tests that pin it

- **`python/tests/test_fanno_exact_composition_tangents_1311.py`:**
  - the grid above, with a bound of 1e-4 (or three times the reference's gap) near the floor and 1e-6 at 300 Pa;
  - linearity in v.
  - Main fails 31 of its 34 cases.
- **`test_fanno_creeping_flow_1216.py::test_adiabatic_static_exit_subfloor_dY_matches_a_difference`:** the former strict xfail of #1303, now a test.
- **`test_fanno_channel_flow_derivatives.py` and `test_fanno_d_ddir_blocked_1203.py`:** rewritten where they asserted d_ddir's stencil behaviour (bit-identity under a round-off component, a 0 off the simplex). They now assert linearity, and the derivative off the simplex equal to the clamped one.
