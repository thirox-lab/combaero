# Exit Dynamic Head Loss at Outflow Boundaries (#1128)

## Problem and Physical Basis

`PressureBoundary.coupling` (`"auto"`, `"static"`, `"total"`) exposes `exit_head_lost(is_outflow)` so that a flow element discharging to an outflow pressure boundary dumps its exit dynamic head into the atmosphere or dump volume (`P_exit = max(P_boundary, P*)`) under `"auto"` and `"static"`, while keeping total-pressure continuity (`Pt_exit = Pt_boundary`) under `"total"` or on inflow (#360, #637).

Before #1128, only `ChannelElement` and `PortDischargeElement` honoured `exit_head_lost`. `PressureLossElement`, `AreaChangeElement`, `LosslessConnectionElement`, and pass-through `OrificeElement` (`bore == D_pipe` or `bore == D_chamber`) always coupled `Pt_exit` to `Pt_boundary`, silently acting as ideal diffusers at outflow boundaries even under `coupling="auto"` or `"static"`. On Probe S4 (`PB 1.05 bar -> 0.5 m x 50 mm channel -> MomentumChamberNode(A) -> X -> PB 1.00 bar`):

| Element `X` (under `coupling="auto"`) | `main` (`e427b904`) `m_dot` [kg/s] | `main` `exit_head_recovered` [Pa] | `#1128` `m_dot` [kg/s] | `#1128` `exit_head_recovered` [Pa] |
|---|---:|---:|---:|---:|
| `ChannelElement` (`L = 1e-4 m`, reference) | 0.19562 | 0.00 | 0.19562 | 0.00 |
| `PressureLossElement` (`zeta = 0`, `area = A`) | 0.42333 | 25567.18 | 0.19562 | 0.00 |
| `LosslessConnectionElement` (from chamber `A`) | 0.42333 | 25567.18 | 0.19562 | 0.00 |
| `AreaChangeElement` (`F0 = A, F1 = A`, pass-through) | 0.42333 | 25567.18 | 0.19562 | 0.00 |
| `OrificeElement` (`bore = D`, pass-through) | 0.42333 | — | 0.19562 | — |
| `AreaChangeElement` (`F0 = A, F1 = 4A`) | 0.24225 | 410.82 | 0.23368 | 0.00 |

`CrossflowSegmentElement`'s station predicate (#1150) is intentionally untouched: a boundary feed remains `entry_K = 0`, and the feed-head family (#1152, #1212, #1153, #1154) is out of scope here.

## Governing Formulation

### 1. Exact Discharge Total Pressure (`solver::discharge_total_pressure`)

Given a boundary static pressure $P_b$, upstream total temperature $T_t$, mass flux $G = |\dot{m}| / A$, and composition $Y$, what upstream stagnation pressure $P_{t,\mathrm{req}}$ is required for a stream of flux $G$ to discharge at static pressure $P_b$ (losing its exit dynamic head) without exceeding Mach 1?

- **Incompressible (`incompressible = true`):**
  $$P_{t,\mathrm{req}} = P_b + \frac{G^2}{2\,\rho(T_t, P_b, X)}$$
  with exact analytic derivatives from ideal-gas $\rho = P_b / (R(X)\,T_t)$:
  $$\frac{\partial P_t}{\partial P_b} = 1 - \frac{q}{P_b}, \quad \frac{\partial P_t}{\partial T_t} = \frac{q}{T_t}, \quad \frac{\partial P_t}{\partial G} = \frac{G}{\rho}, \quad \frac{\partial P_t}{\partial Y_k} = \frac{q}{R}\frac{\partial R}{\partial Y_k}$$
  where $q = G^2 / (2\rho)$ and $\partial R / \partial Y_k$ is computed analytically via `specific_gas_constant_dY` (zero finite differences).

- **Compressible (`incompressible = false`):**
  Let $(P^*, P_t^*)$ be the sonic static and total pressure of mass flux $G$ at $(T_t, X)$ (`solver::sonic_static_pressure` + isentropic relation at $M = 1$).
  - **Subsonic ($P_b > P^*$):** $P_{t,\mathrm{req}}$ is the unique subsonic stagnation pressure whose isentropic static state at $(P_b, T_t, G, X)$ has Mach $M < 1$, obtained directly from `solver::momentum_chamber_residual_and_jacobian(P_b, P_b, m = G, T_t, Y, A = 1)` (the exact real-gas $c_p(T)$ isentropic relation already used by `MomentumChamberNode`). All four partials ($\partial P_t / \partial P_b$, $\partial P_t / \partial T_t$, $\partial P_t / \partial G$, $\partial P_t / \partial Y_k$) come directly from the exact analytic derivatives of that closure.
  - **Choked ($P_b \le P^*$):** The exit chokes at $M = 1$ with under-expanded exit static pressure $P^* \ge P_b$, so
    $$P_{t,\mathrm{req}} = P_t^*(T_t, G, Y) = P^*(T_t, G, Y)\,\exp\!\left(\Delta s_{T^* \to T_t} / R\right),$$
    which is independent of $P_b$ ($\partial P_t / \partial P_b = 0$).
  - **$C^1$ continuity across the sonic choke (envelope theorem):** At fixed $(T_t, G, Y)$, $P_t(P, T_t, G, Y)$ attains its strict minimum at the sonic pressure $P = P^*$ where $\partial P_t / \partial P \big|_{P=P^*} = 0$. Therefore:
    $$\frac{\partial P_t^*}{\partial \theta} = \left.\frac{\partial P_t}{\partial \theta}\right|_{P=P^*} + \underbrace{\left.\frac{\partial P_t}{\partial P}\right|_{P=P^*}}_{=\,0}\frac{\partial P^*}{\partial \theta} = \left.\frac{\partial P_t}{\partial \theta}\right|_{P=P^*}$$
    for $\theta \in \{T_t, G, Y_k\}$. Evaluating `momentum_chamber_residual_and_jacobian` at $P = P^*$ yields the exact analytic sonic partials and guarantees $C^1$ smoothness of $P_{t,\mathrm{req}}$ and all derivatives across $P_b = P^*$ with zero smoothing constants and zero finite differences.

### 2. Element Residuals

- **`PressureLossElement`:** In each flow direction (`_directed`), when the flow-downstream node declares `exit_head_lost(True)` and the element has a velocity area (`_exit_head_area() is not None`, i.e., explicit/inherited area or correlation area, excluding the `DEFAULT_CHAMBER_AREA` placeholder), the residual becomes $P_{t,\mathrm{req}}(P_{b,\mathrm{dn}}, T_{t,\mathrm{up}}, |\dot{m}|/A, Y_{\mathrm{up}}) - P_{t,\mathrm{up}}(1 - \xi) = 0$. Forward and reverse legs blend smoothly through zero flow via `_flow_blend`.
- **`LosslessConnectionElement`:** When the flow-downstream node declares `exit_head_lost(True)`, uses `self.area` (compressible) if set, else `up_node.area` (and `up_node.regime`) if leaving a `MomentumChamberNode` (Option A decision on #1128). Leaving a `PlenumNode` with no area keeps `Pt_in - Pt_out = 0` bit-identically.
- **Pass-through `OrificeElement` (`bore == D_pipe` or `bore == D_chamber`, including plates):** Shares `_discharge_total_row` and mirrors `LosslessConnectionElement._exit_head_spec` (using `up_node.area` and `up_node.regime` when leaving a `MomentumChamberNode`, else `self.area` and `self.regime`), preserving #1016's exact equivalence with `LosslessConnectionElement`.
- **`AreaChangeElement`:**
  - Normal (`F0 != F1`): C++ `area_change_flow` and `area_change_element_flow` gain `exit_head_lost` (`exit_head_lost_b` / `exit_head_lost_a`), adding the exact exit dynamic head ($P_t(M_e) / P_e$ in compressible flow, $q_e$ in incompressible flow) in the Brent target equation using `DualN` forward-mode automatic differentiation for all exact partials. When low-Re viscous loss causes the exit face to reach $M_e = 1$ before the inlet reaches $M_1 = 1$, the upper bracket clamps cleanly to the exit-face sonic Mach $M_{1,\max}$ (with exact `DualN` derivative of the $M_e = 1$ implicit constraint).
  - Pass-through (`F0 == F1`): Uses `_discharge_total_row` with `(F0, regime == "incompressible")`.

### 3. Design Alternatives Considered (Option A vs Option B)

Two interpretations of `MomentumChamberNode -> LosslessConnectionElement(area=None) -> PressureBoundary(coupling="auto")` were evaluated and put to the owner on issue #1128:

- **Option A (chosen by owner):** `PressureBoundary(coupling="auto")` at an outflow always means the static pressure there, whatever link connects it (`PressureLossElement` with area, `AreaChangeElement`, `LosslessConnectionElement(area=None)` leaving a `MomentumChamberNode`, and pass-through `OrificeElement`). Leaving a `PlenumNode` (zero velocity) is unchanged. Harnesses that intentionally impose stagnation pressure `Pt` at a port switch that boundary to `coupling="total"`.
- **Option B (rejected):** Leave `LosslessConnectionElement(area=None)` as pure `Pt_in = Pt_out` even out of a `MomentumChamberNode`. Rejected because it leaves the #1128 inconsistency unresolved on `chamber -> Lossless -> PB(auto)` and forces users to know about `PortDischargeElement` to avoid an implicit ideal diffuser.

### 4. Solver Robustness Findings (Root-Caused, Not Masked)

1. **`test_fraction_unanchored_1004` (`NO_PROGRESS` stall):** Root-caused to a pre-existing `hybr` failure on blocked `ConstantFractionLoss` / `LinearThetaFractionLoss` networks without `m_ref` (4/35 `(zeta, P_a)` cases already failed identically on `main` `e427b904`). Outside `FLOW_REVERSAL_M_SCALE = 1e-4 kg/s`, `dF/dm = 0`, so `hybr` overshoots across `m = 0` onto the opposite branch cliff and aborts after 10 rejected trust-region steps (`NO_PROGRESS`), while `NetworkSolver` skipped its existing LM fallback because `PressureLossElement` had no `regime="compressible"` attribute. Gating ONLY the post-failure LM fallback (never `hybr`'s 65% budget cut or stall handover, and never inside `_in_warmstart_proxy`) on `_has_unanchored_fraction` (`solver.py`) rescues 70/70 swept cases and 10/10 `p.Pt` seeds to `|F| ~ 1e-11`.
2. **Static-discharging fraction-loss cold seed (`solver.py`, `random` seed 35):** In `_resistance_network_guess`, a fraction loss's Picard ratio row previously enforced flow-independent `Pt_dn - (1 - xi) Pt_up = 0` (`L[col, col] = 0`) even when the element discharged with an exit area into a static boundary (seeding `random` 35 at 34.5 kg/s vs its 0.238 kg/s choked root). Carrying the exit dynamic-head linearization `mode[k] * q_dir * (|m[k]| + SEED_MDOT_FLOOR)` (`q_dir = 1 / (2 rho A^2)`) on `L[col, col]` preserves `#989`'s fraction-ratio law while respecting exit head (`random` 35 solves cold in `0.31 s`, `random` 4 in `1.48 s`).
3. **Legacy-tee forward-start retry on failed reversed iterates (`solver.py`, `legacy_tee` seed 34):** Under static discharge, a merging `TeeJunctionElement` cold solve that strays into reversed flow can stall past sonic (`no_progress`) instead of converging to a refused reversed root (`inconsistent`). `_retry_legacy_tees_forward` now scans all `TeeJunctionElement`s at the failed solve's best iterate (not only `__inconsistent_elements__`) and also runs before `homotopy_retry` / `bc_ramp_retry` so those retries cannot burn the budget first, restoring `legacy_tee` 34's forward root in `0.38 s`.
4. **Exact-root guard (`solver.py`):** In `solver.py`, MINPACK's progress test can stop at machine precision when the raw Pa norm is `~1e-11` (above `1e-12`) even though the scaled dimensionless norm is `< 1e-12` (`~1e-16`). Accepting `best_scaled_norm <= _EXACT_ROOT_TOL and final_norm <= _RESIDUAL_TOL` closes this exact-root false negative without loosening `_RESIDUAL_TOL`.

## Verification and Test Pointers

- C++ kernel tests (analytic vs FD, sonic $C^1$ handover, low-Re exit-face choke): [tests/test_exit_head_lost_1128.cpp](../../tests/test_exit_head_lost_1128.cpp).
- Python regression, Probe S4, pass-through orifice equivalence, multi-seed blocked-fraction, reversal, and full-network FD Jacobian tests: [python/tests/test_exit_head_lost_1128.py](../../python/tests/test_exit_head_lost_1128.py), [python/tests/test_exit_head_recovered.py](../../python/tests/test_exit_head_recovered.py), [python/tests/test_orifice_bore_chamber_1016.py](../../python/tests/test_orifice_bore_chamber_1016.py).
- Cross-source Behrendt (2008) pins (re-pinned with unchanged band widths after chamber static sits at `p_h` rather than `p_h - q_gas`): [python/tests/test_effusion_network_runner.py](../../python/tests/test_effusion_network_runner.py).
