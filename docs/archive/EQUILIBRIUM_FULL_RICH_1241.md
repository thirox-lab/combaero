# Full equilibrium of rich products (#1241): provenance record

0.13.0, theme 1 of the release plan (#1341). Owner decisions: "methanation
always" (2026-11-24; written as "1243", meant #1241), NH3 synthesis always
(plan, 2026-11-25); the representable bound, and the per-stream one-sided
slope, on 2026-11-27.

## Formulation

- **Which reactions run.** `reforming_species_present` in `src/equilibrium.cpp`.
  - Reforming, with methanation as its backward direction, runs where the burned state holds a reducing species (H2, CO, a hydrocarbon or NH3) and both carbon and hydrogen.
  - NH3 synthesis and decomposition run where it holds one and both nitrogen and hydrogen.
  - There is no threshold, so nothing switches as a trace goes to zero.
  - Without a reducing species neither reaction can move, and the shift-only path is kept as it was: lean and fuel-free states stay bit-identical, value and Jacobian.
- **The representable bound.** `joint_key_unrepresentable`, `drop_unrepresentable_keys`.
  - A reaction's key (the species no other reaction holds) whose equilibrium amount, with every other amount held, is below the smallest normal double is held at exactly 0. Its reaction is run to completion away from it.
  - This is applied only after a failed solve, so a state that converges is unchanged. Seeding skips such a key.
  - In the Jacobian, a fed amount of such a key goes through its reaction completely. The burned amounts become R n0, with R's key column e_k - nu/nu_k.
- **One-sided slopes at zero flow.** `dX_dir` of `combustion_equilibrium_and_jacobian`; the zero-flow override in `equilibrium_from_streams` in `src/solver_interface.cpp`.
  - Where the burned state lacks a species an arriving reactant feeds, the reaction set and the trace ladder are taken along that reactant's direction.
  - Where the direction feeds a species no reaction of the solved basis holds, the plain reaction set is taken instead (carbon arriving in carbon-free air or in rich H2).
  - The combustor passes each zero-flow stream's mixing direction. Every other column, and every stream with flow, keeps the matrix.
- **The limit.** `REFORMING_JACOBIAN_LIMIT_RTOL`.
  - From the first converging trace, the ladder steps down by 100 until two traces' derivatives agree to 1e-8: along the direction if one is given, else over the matrix, each relative to its own infinity-norm.
  - The kept derivatives are the larger trace's, so where the first trace is already in the limit's regime J is as before.
  - The step-down runs only along `dX_dir`. Without a direction, the seed is the same trace of every fed species, which is no reactant's direction, so its limit means nothing.
  - The floor is sqrt(DBL_MIN): a species quadratic in the trace (N2 beside the NH3 it forms) stays normal.
- **Light-first fallback.** After the bound, a solve that still fails:
  - every hydrocarbon is reformed completely with steam;
  - the shift and NH3 are solved without reforming, by the full solver;
  - each hydrocarbon is added back at its own reforming equilibrium at that converged state, or left at 0 below the normal range;
  - the set is solved once more, without the fallback.

  Representability is judged at a converged state, not at the stalled iterate. The nightly search found the case: CH4 and heptane smoothing traces of 8e-18 at 500 K and 0.1 bar.
- **The NH3 start.** `solve_nh3_decomposition_isothermal` returns its bracketed iterate at the iteration cap instead of throwing. It is a start; the joint solve decides convergence.

## Ground truth

- **Values.** Cantera 3.2 HP on combaero's own species set (NH3 included) and NASA-9 data. Rich states where Cantera's O2 stays below 1e-7 (155 random draws): worst |dT| 9.7e-4 K, against 5.43 K before. Pinned cases are in `python/tests/test_equilibrium_full_rich_1241.py`.
- **Derivatives.** Forward differences along the direction, where a difference resolves the limit (2e-5), and continuity, J(0) = J(1e-40 kg/s of air) on an ordinary state without the ladder (1e-7).
- **The NH3 trace.** ln Q - ln K of the decomposition, from the public h and s of the pure species (1e-14).

## Why not the alternatives

- **The trace-limit route of 725e4f57.** It reformed a hydrocarbon trace completely and then applied the shift, with a 1e-6 margin constant.
  - Its "frozen" route kept the shift alone wherever the burn holds CO or H2 of its own, which is "methanation always" reversed.
  - It gated on a fed hydrocarbon, so an H2 trace in air had no route.
  - The representable bound reaches the same limit inside the joint solve, with no second route and no margin constant.
- **Keeping a cut for traces.** Rejected by the owner: the cut is #1241's step.
- **A finite-difference slope at zero flow.** Not allowed: no hidden finite difference in a kernel.
- **A Jacobian matrix at zero flow.** None is right in dry air. The slopes along CH4 and CO do not add, by 11 % at 335 K and typically 0.2-5 %, because steam is a trace too. In humid air they add exactly.

## Choices not forced by the spec

- **The bound at DBL_MIN (normal range), not denorm_min.** With denorm_min, n-hexane at 1.4e-320 (four digits) quantised its row's residual to 1e-5 and stalled the propane-to-pentane rows. The NH3 start's cap is a separate question.
- **The gap metric relative to each vector's infinity-norm.** A trace species far below the largest entry can converge slowly unseen, by design: its effect on T and the majors is that small.
- **Step factor 100 and stop on agreement.** The same factor as the ladder's way up.

## Dead ends

- **The union predicate, every column of dXb/dX counting as fed.** It did not fix the fuel-free columns, which need a direction, and it refused two NH3 burns near phi 1.
- **The fixed first trace of 1e-20.** It lies past the bend for fuel-only blends: dY_NH3/dm_air was 75.5 where the limit is 1594.

## Open after this change

- **#1343:** 5 of 4000 extreme states are refused (near-pure CO with a hydrogen trace and about 1e-40 of air; CH4 at 63 K).
- **#1347:** with an NH3 trace (1e-17), dT/dX in a cold rich isobutane burn is noisy at 1e-4. Without N2 it is exact. One nightly cell sits at 1.01e-4 against its 1e-4 bound.
- **#1292, #1345:** OPEN sweep patterns, the one-sided limits at zero fuel and zero air.

## What a difference cannot check

The one-sided limit at zero flow can lie far below any step a solver takes. For air reaching a fuel-only blend at 300 K, the N2 goes completely to NH3, and the first oxygen's products reform beside the fuel, only below about 1e-26 kg/s of air. There J(0) is the partial-oxidation slope, 1.53e4 K per kg/s, while a difference over 1e-15 to 1e-11 kg/s sees 1.39e6. Cold fuel-free states (below 500 K) bend the same way, the #1292 pattern. The owner accepted this compromise for edge cases: the limit is the defined, checkable quantity.

## Measured

- **Bit identity against main e427b904** (800 draws, seed 20261127): lean 200/200 and fuel-free 200/200 identical in value and Jacobian.
- **Wide search** (4000 draws, seed 5, T_in 60-2000 K, air 1e-40 to 1): 74 states main refuses now converge; 5 converge on main only. Those 5 are near-pure CO with a hydrogen trace and about 1e-40 of air, or CH4 at 63 K; filed as #1343.
- **Warm fuel traces in air** (1e-20 to 1e-13, 300-900 K): 0 refused, against 151 before the bound (main froze them below its cut).
- **#1290:** converges.

## Tests that pin it

- `test_equilibrium_full_rich_1241.py`: the Cantera pins, the warm-trace grid, the directional slopes, and that dX_dir is a no-op elsewhere.
- `test_combustion_jacobians.py`:
  - `test_a_fuel_only_equilibrium_combustor_moves_with_the_air_flow` holds the continuity check, which pins the limit tolerance: at 1e-5 the gap is 3.3e-6.
  - `test_heptane_with_air_below_the_steam_floor_converges` covers #1290.
  - `test_an_unconverged_equilibrium_is_refused` now uses #1343's case.
- `test_equilibrium_nh3_1314.py::test_the_lean_softmin_trace_of_nh3_is_decomposed_to_equilibrium`.
- **Falsified:**
  - with the bound off, 92 of 110 guard tests fail;
  - on main, the Cantera pins, the directional tests and 4 of 8 #1339 cases fail;
  - with the zero-flow override off, all 32 zero-air combustor tests fail;
  - with the plain-set fallback off, 10 of 16 directional tests fail;
  - with the R transform off, 1 fails (heptane arriving in carbon-free air).
- **The sweep** (`_jacobian_sweep.py`; families combustor, walled_combustor, fuel_controller, film):
  - main had 8 mismatches, all at exact zeros; 7 of those 8 are gone;
  - this branch has 15, all at exact zeros: zero fuel when cold (#1292) and zero air (#1345), the one-sided limits the ladder cannot reach;
  - both are OPEN patterns, capped at about twice the measured error.
- **#1339** (cold rich isobutane with an O2 trace): fixed. The branch J took there was the cut.
- **C++:** four reforming tests in `test_thermo_transport.cpp` used N2 as the inert reference. They now normalise by N2 + NH3/2.
