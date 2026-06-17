# Systems

This folder holds the **system builders** — the functions that construct the
PowerSystems.jl (`PSY.System`) objects that every model, simulation, and example
script in this repository consumes. Each subfolder targets one test system and
returns a `System` that is ready to feed into a UC / ED problem template.

Three systems are provided, in increasing order of size and complexity:

| System | Folder | Buses (approx.) | Use it for |
| ------ | ------ | --------------- | ---------- |
| 5-bus  | `5bus/` | 5   | Fast iteration, transmission-expansion / flow-cancelling experiments |
| RTS    | `RTS/`  | 73  | Realistic-but-small benchmark; circular-flow and loss studies |
| CATS   | `CATS/` | ~8500 | Large California system; production-scale stress tests |

All three integrate with the script libraries in
[`../SiennaScripts`](../SiennaScripts/README.md).

---

## `5bus/build_5bus.jl`

Matpower 5-bus base system plus tooling for transmission-expansion / candidate-line
studies (used by the flow-cancelling work).

Key entry points:

| Function | Purpose |
| -------- | ------- |
| `build_matpower_5bus_with_updated_lines()` | Build the 5-bus system, attach single-time-series load, linearize generation costs, and adjust line ratings. The standard starting point. |
| `update_generation_costs!(sys)` | Replace thermal cost curves with linear (no fixed / no-load) costs. |
| `update_line_ratings!(sys, thermal_multiplier)` | Scale line ratings and pin a couple of specific lines. |
| `add_candidate_line_data!(sys)` / `add_candidate_line_data_without_parallel!(sys)` | Add candidate (investment) lines, with or without modeling them as parallel circuits. |
| `add_new_line_without_parallel!`, `update_line_to_parallel!`, `candidate_projects_data`, `reconductoring_candidate_data` | Lower-level helpers for building candidate-line and reconductoring scenarios. |

## `RTS/build_rts.jl`

Modified RTS-GMLC day-ahead system. This is the main system used by the
**CircularFlows** examples.

Key entry points:

| Function | Purpose |
| -------- | ------- |
| `build_rts_system()` | Build the RTS system, transform to hourly single time series, and configure the HVDC line (lossless, ±limits). |
| `set_renewable_costs!(sys, cost_sign)` | Set renewable cost curves with a chosen sign. A **positive** sign makes the HVDC loop profitable (circular flow appears); a **negative** sign removes the incentive. This sign trick drives the Scenario A / B contrast in the examples. |
| `make_uc_template()` | PTDF-based UC `ProblemTemplate` with an AC power flow evaluation (loss + voltage-stability factors). |
| `build_rts_model_with_quadratic_losses(sys)` | RTS model with the quadratic `I²R` loss term added to the copper-plate balance. |
| `build_rts_uc_models_hv(...)` / `build_rts_ed_models_hv(...)` | Device models with a `filter_function` that removes low-voltage branches from the thermal-limit constraints (PTDF still full-size for losses). |
| `build_rts_ed_models_acopf()` | ED device models for the `ACPPowerModel` (ACOPF) formulation. |

## `CATS/build_cats.jl`

Large California Test System (CATS). The system itself is stored as serialized
artifacts in this folder and loaded from disk rather than rebuilt.

Serialized data (large binaries):

- `CATS_saved_sys.json` / `CATS_saved_sys_metadata.json` / `CATS_saved_sys_time_series_storage.h5` — full system.
- `CATS_saved_reduced_sys.*` — reduced variant.

Key entry points:

| Function | Purpose |
| -------- | ------- |
| `build_cats_system(cats_json_path)` | Deserialize a CATS `System` from the JSON path. |
| `set_cats_renewable_costs!(sys, cost_sign)` | Same renewable-cost sign trick as RTS, for inducing/suppressing circular flows. |
| `scale_cats_loads!(sys, scale_factor)` | Uniformly scale loads. |
| `build_cats_uc_models_hv(...)` / `build_cats_ed_models_hv(...)` | CATS device models with the low-voltage branch `filter_function`. Note CATS uses `Transformer2W` (not `TapTransformer`) and has `SynchronousCondenser`. |
| `build_cats_ed_models_acopf(...)` | CATS ED device models for the ACOPF formulation. |
| `CATS_UC_MODELS` / `CATS_ED_MODELS` | Default device-model dictionaries for CATS. |

> The HV `filter_function` removes the per-branch flow-limit constraints for
> low-voltage branches but keeps their contribution to the quadratic loss term:
> the full system PTDF must still be passed to the PSI `NetworkModel`.
