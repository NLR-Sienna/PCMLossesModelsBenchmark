# CircularFlows

Tools to **detect circular (loop) flows** in a solved power system. After an
optimization (UC, ED, or ACOPF) is solved, the branch flows are assembled into a
directed weighted graph; any directed cycle in that graph is a circular flow — power
chasing itself around a loop, often induced by an HVDC link plus a cost incentive.

The recurring demonstration uses two scenarios:

- **Scenario A** — renewable cost > 0 (or HVDC loop enabled) ⇒ the loop is profitable
  ⇒ a circular flow is detected.
- **Scenario B** — renewable cost < 0 (or HVDC disabled) ⇒ no incentive ⇒ no circular
  flow.

This folder builds on [`../../Systems`](../../Systems/README.md) (RTS / CATS builders)
and the model/simulation library in [`..`](../README.md).

---

## Core module files

| File | Purpose |
| ---- | ------- |
| `circular_flows.jl` | The detection engine. Builds a flow graph and finds cycles. Multiple graph sources: `build_graph_from_ptdf`, `build_graph_from_pf_aux_variables`, `build_graph_from_acopf_variables`; HVDC edges via `add_hvdc_edges!`; cycle search via `find_circular_flows`. Also `read_hvdc_flow_variables` and `build_graph`. |
| `mapped_indices.jl` | `MappedIndices <: AbstractArray` — maps internal graph node indices (from `simplecycles`) to PSY bus numbers in one indexing step. Helpers: `lookup_arrays`, `e2i`, `i2e`. |
| `print_utils.jl` | Reporting. `compare_voltage_stability_factors` / `print_stability_comparison` build per-bus comparison tables between two scenarios (with vs without circular flows). |
| `rts_example.jl` | End-to-end RTS validation (lossless PTDF UC, HiGHS), Kestrel-style activation. |
| `rts_example_quadratic_losses.jl` | End-to-end RTS build validation with the quadratic loss term. Requires Gurobi with `NonConvex=2` to fully solve. |

---

## `scripts/` catalog

Each script is a self-contained example. The `_local` variants activate the repo
environment and are meant to be run **from the repo root** (e.g.
`julia SiennaScripts/CircularFlows/scripts/rts_example_local.jl`); the `_kestrel`
variants are deployed at the project root on the Kestrel HPC cluster.

| Script | System | Stages / formulation | Target | Solver |
| ------ | ------ | -------------------- | ------ | ------ |
| `rts_example_local.jl` | RTS | Lossless PTDF UC | local | HiGHS |
| `rts_example_kestrel.jl` | RTS | Lossless PTDF UC | Kestrel | HiGHS |
| `rts_example_quadratic_losses_local.jl` | RTS | UC + quadratic `I²R` losses | local | Gurobi (`NonConvex=2`) |
| `rts_example_quadratic_losses_kestrel.jl` | RTS | UC + quadratic `I²R` losses | Kestrel | Gurobi (`NonConvex=2`) |
| `rts_example_sim_quadratic_losses_local.jl` | RTS | UC (PTDF) → ED (quadratic losses) via `SemiContinuousFeedforward` | local | HiGHS + Ipopt |
| `rts_example_sim_filter_quadratic_losses_local.jl` | RTS | Same as above, with low-voltage branch filtering | local | HiGHS + Ipopt |
| `rts_example_sim_acopf_ed_local.jl` | RTS | UC (PTDF + HV filter) → ED (`ACPPowerModel`, all lines) | local | HiGHS + Ipopt |
| `cats_example_local.jl` | CATS | Lossless PTDF UC | local | HiGHS |
| `cats_example_kestrel.jl` | CATS | Lossless PTDF UC | Kestrel | HiGHS |
| `cats_example_quadratic_losses_local.jl` | CATS | UC (PTDF) → ED (quadratic losses); binaries fixed by feedforward ⇒ pure NLP | local | HiGHS + Ipopt |
| `cats_example_quadratic_losses_kestrel.jl` | CATS | UC (PTDF) → ED (quadratic losses) | Kestrel | HiGHS + Ipopt |
| `cats_example_sim_quadratic_losses_local.jl` | CATS | UC → ED (quadratic losses) via `SemiContinuousFeedforward` (no manual binary fixing) | local | HiGHS + Ipopt |
| `cats_example_sim_filter_quadratic_losses_local.jl` | CATS | Same, with low-voltage branch filtering to reduce model size | local | HiGHS + Ipopt |
| `cats_example_sim_acopf_ed_local.jl` | CATS | UC (PTDF + HV filter) → ED (`ACPPowerModel`, all lines) | local | HiGHS + Ipopt |
| `test_rts_filter_quadratic_losses.jl` | RTS | Validation test for the filtered-line quadratic-loss simulation (asserts Scenario A finds ≥1 loop, B finds 0) | local | HiGHS + Ipopt |

### Notes on the variants

- **`sim_*`** scripts use a PSI `Simulation` with `SemiContinuousFeedforward`, so UC
  commitment binaries are propagated to ED automatically — ED becomes a pure NLP and
  no MIQP solver (Gurobi) is needed.
- **`filter`** scripts exclude low-voltage branches from the thermal-limit
  constraints (`NetworkFlowConstraint`) to shrink the model, while still using the
  **full** system PTDF for the loss term — losses on filtered LV branches are still
  captured.
- **`acopf_ed`** scripts solve ED as an `ACPPowerModel`. They expose two modes via
  `ignore_pf_ed`:
  - `true` (default) — circular flows read directly from ACOPF optimization variables.
  - `false` — a post-solve AC power flow runs; circular flows come from PF aux
    variables and voltage stability factors are available from PSI results.

---

## Start here

For a first run with no special solver licensing, use:

```bash
julia SiennaScripts/CircularFlows/scripts/rts_example_local.jl
```

It builds the RTS system, solves a lossless PTDF UC with HiGHS, and prints the
detected circular flow for Scenario A (and its absence for Scenario B). From there:

- Add losses → `rts_example_sim_quadratic_losses_local.jl` (still HiGHS + Ipopt).
- Try ACOPF ED → `rts_example_sim_acopf_ed_local.jl`.
- Scale up to CATS → the `cats_example_*_local.jl` family.

**Solver requirements at a glance:** HiGHS (open source) for lossless and
feedforward simulations; Ipopt for the NLP ED stage; Gurobi with `NonConvex=2` only
for the non-simulation quadratic-loss MIQP scripts (`*_quadratic_losses_local/kestrel.jl`).

## Background

Design notes for this module live in
[`../../docs/superpowers/plans/`](../../docs/superpowers/plans/) — see the
`circular-flows-module`, `circular-flows-unify-api-and-examples`,
`acopf-ed-circular-flows`, and `cats-filtered-lines-quadratic-losses` plans.
