# SiennaScripts

This is the **reusable model and simulation library** for the benchmark. The files
at this level are not meant to be run directly — they define the functions that the
example scripts (in the subfolders and at the repo root) call to build, solve, and
post-process Sienna optimization problems.

## Layered design

```
Systems/            build the PSY.System (network + devices + time series)
   │
SiennaScripts/      build & solve DecisionModels / Simulations, extract results
   │
   ├── CircularFlows/   detect loop flows in a solved system   (see its README)
   └── FlowCancelling/  transmission-investment + flow cancelling (see its README)
```

Most workflows: pick a builder from [`../Systems`](../Systems/README.md), pass the
`System` to a `build_*` / `run_*` function here, then read results with helpers from
`utils.jl`.

---

## Major files

### `build_models.jl` — single `DecisionModel` construction
Builds one optimization model (typically UC) and the loss-augmentation machinery.

| Function | Purpose |
| -------- | ------- |
| `make_ptdf_model_without_losses(...)` / `build_ptdf_model_without_losses(...)` | Construct (and build) a lossless PTDF model. |
| `make_acopf_model(...)` | Construct an ACOPF (`ACPPowerModel`) model. |
| `add_current_loss_variables!`, `add_current_loss_to_copperplate_balance!`, `add_quadratic_current_loss_to_copperplate_balance!` | Add linear / quadratic `I²R` loss terms to an existing model. |
| `add_current_loss_constraint_approximation!`, `add_current_loss_constraint_quadratic_approximation!` | Iterative loss-constraint approximations around a previous solution. |

Also defines `DEFAULT_UC_MODELS`, `DEFAULT_ED_MODELS`, and the optimizer constants
`DEFAULT_MILP_OPTIMIZER` (Xpress) and `DEFAULT_NLP_OPTIMIZER` (Ipopt), which fall
back to `nothing` when the solver package is not loaded.

### `build_simulations.jl` — multi-stage UC→ED `Simulation` construction
Each function returns a configured (unexecuted) `Simulation`.

| Function | Purpose |
| -------- | ------- |
| `build_uc_ed_simulation_with_no_losses(...)` | Baseline lossless UC→ED. |
| `build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses(...)` | Linear losses in UC, quadratic in ED. |
| `build_uc_ed_simulation_with_acopf(...)` | ED solved as ACOPF. |
| `build_uc_ed_simulation_with_acopf_and_uc_linear_losses(...)` | ACOPF ED + linear-loss UC. |
| `build_uc_double_ed_simulation_with_acopf(...)` / `..._and_uc_linear_ed_quadratic_losses(...)` | Two ED stages (e.g. forecast + AC re-dispatch). |
| `build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(...)` | Quadratic ED losses without voltage-dependent terms. |

Two flags recur: `ignore_pf_uc` / `ignore_pf_ed` skip the post-solve AC power flow
at each stage (recommended `true` for UC; `false` for ED keeps PF aux variables such
as loss factors and voltage stability factors available).

### `run_models.jl` — build + solve + extract (single model)
| Function | Purpose |
| -------- | ------- |
| `run_lossless_model(sys)` | Build + solve the lossless PTDF UC model; return model, results, bus injections. |
| `run_linear_loss_model(...)` | One linear-loss pass. |
| `run_iterative_linear_loss_model(...)` | Iterate the linear-loss approximation to convergence. |

### `run_simulations.jl` — build + execute + extract (simulations)
| Function | Purpose |
| -------- | ------- |
| `run_uc_ed_lossless_simulation(...)` | Execute the lossless UC→ED simulation. |
| `run_uc_ed_acopf_simulation(...)` | ACOPF ED. |
| `run_uc_ed_quadratic_loss_simulation(...)` / `run_iterative_uc_ed_quadratic_loss_simulation(...)` | Quadratic-loss ED, single pass / iterated. |
| `run_uc_linear_loss_ed_acopf_simulation(...)` / `run_iterative_uc_linear_ed_acopf_simulation(...)` | Linear-loss UC + ACOPF ED, single / iterated. |

### `utils.jl` — result post-processing
Helpers that read PSI results and return plain vectors/matrices (axis-1 = component,
axis-2 = timestep), avoiding DataFrame juggling: `get_bus_ax`, `get_bus_injection`,
`get_bus_loss_factors`, `get_total_R_per_arc`, `get_total_X_per_arc`,
`get_power_flow_voltage_mag`, `get_power_flow_arc_voltage_mag`, `get_RX_vector`.

### `add_hvdc.jl`
`add_internal_hvdc!(sys)` injects two synthetic `TwoTerminalGenericHVDCLine`
components into the CATS network, creating a closed loop that can carry circular
flows — the setup used by the CircularFlows CATS examples.

### `formatter/`
JuliaFormatter configuration for the repository.

---

## Subfolders

- **[`CircularFlows/`](CircularFlows/README.md)** — detect loop (circular) flows in
  solved RTS / CATS systems, with a large catalog of runnable examples.
- **[`FlowCancelling/`](FlowCancelling/README.md)** — transmission-investment models
  with Big-M flow-cancelling correction terms.
