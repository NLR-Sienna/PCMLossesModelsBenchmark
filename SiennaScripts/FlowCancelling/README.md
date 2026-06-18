# FlowCancelling

Models for **transmission investment with flow-cancelling correction terms**. The
formulation adds candidate (investment) lines to a PTDF network model and uses Big-M
linking so that a candidate line only shifts power flows when it is built. Optional
quadratic `I²R` losses can be layered on top.

The idea: a candidate line's effect on existing-line flows is represented by a PTDF
**shift term** (`BranchCancellingFlowVariable`). A Big-M constraint forces that shift
to zero unless the line's investment binary (`BranchInvestmentVariable`) is on — so an
un-built line "cancels" its own flow contribution.

This folder pairs with the 5-bus candidate-line tooling in
[`../../Systems/5bus/build_5bus.jl`](../../Systems/README.md).

---

## `build_models.jl`

Custom PSI types:

| Type | Role |
| ---- | ---- |
| `GenerationInvestmentVariable` | Generation investment decision. |
| `BranchInvestmentVariable` | Candidate-line build decision (binary). |
| `BranchCancellingFlowVariable` | Per-candidate PTDF flow shift/cancelling term. |
| `BigMConstraint` | Links the cancelling flow to the build decision. |
| `LineLossTotalApproximation`, `LineLossConstraintApproximation`, `PTDFBranchFlowWithFC` | Quadratic-loss variable / constraint / augmented-flow expression. |

Top-level builders:

| Function | Purpose |
| -------- | ------- |
| `make_base_ptdf_model(sys; ...)` / `build_base_ptdf_model(sys; ...)` | Create (and build) the base PTDF model. `ignore_pf = false` adds an AC power flow with loss-factor calculation. |
| `build_model_with_flow_canceling_terms(sys; ignore_pf = true, device_models = DEFAULT_FC_DEVICE_MODELS, optimizer = DEFAULT_MILP_OPTIMIZER)` | Full flow-cancelling investment model. |
| `build_model_with_flow_canceling_and_quadratic_losses(sys; device_models = DEFAULT_FC_DEVICE_MODELS)` | Flow-cancelling model plus quadratic loss constraints. |

Pass `device_models = RTS_FC_DEVICE_MODELS` for the RTS system (it adds `TapTransformer`,
renewables, hydro and HVDC). Branch types must use `StaticBranch` so the `FlowRateConstraint`
containers exist.

Supporting builders (called by the top-level ones): `add_branch_investment_variables!`,
`add_branch_cancelling_flow_variables!`, `add_bigM_linking_constraints!`,
`add_shift_terms_to_existing_line_constraints!`,
`add_shift_terms_to_candidate_line_constraints!`,
`add_candidate_line_investment_costs!`,
`add_candidate_generation_investment_constraints!`, and the quadratic-loss helpers
`_fc_add_loss_variables!`, `_fc_add_loss_to_copperplate_balance!`,
`_fc_add_quadratic_loss_constraints!`.

### Parallel lines (`-double_circuit`) and arc resolution

On networks with parallel branches (e.g. RTS), `PTDFPowerModel` aggregates parallel circuits
into a single reduced branch named `<name>-double_circuit`. These reduced names — not the PSY
component names — key the `FlowRateConstraint` / `PTDFBranchFlow` axes. The flow-cancelling
code therefore resolves every monitored-branch arc through
`ptdf.network_reduction_data.name_to_arc_map` (`_branch_arc_map` / `_resolve_branch_arc`) and
iterates **constraint axes** rather than PSY components, so each parallel pair is handled once.
`_get_ptdf_shift_term` takes the resolved `(from,to)` arc tuple directly. Candidate lines are
added via the intermediate-bus trick and are never parallel.

> Requires a MILP solver — `DEFAULT_MILP_OPTIMIZER` is Xpress (guarded so the file still loads
> without Xpress; pass `optimizer = ...` to use another MILP solver such as HiGHS).

---

## How to run

There is no `scripts/` subfolder here; the runnable examples live at the repository
root and drive the 5-bus system:

| Script (repo root) | Purpose |
| ------------------ | ------- |
| `script_5bus_flowcancelling.jl` | Minimal flow-cancelling investment run on the 5-bus system. |
| `script_5bus_flowcancelling_yc.jl` | Variant of the above. |
| `tutorial_losses_and_flow_cancelling.jl` | Walkthrough combining losses with the flow-cancelling formulation. |
| `script_rts_flowcancelling.jl` | Flow-cancelling on RTS GMLC (parallel `-double_circuit` lines + TapTransformer); base + stressed scenarios with sanity-check assertions. |

Typical 5-bus flow: build the system with candidate lines
(`build_matpower_5bus_with_updated_lines` + `add_candidate_line_data!` from
`Systems/5bus/build_5bus.jl`), then pass it to `build_model_with_flow_canceling_terms`.

Typical RTS flow: `build_rts_system()` + `add_candidate_lines_rts!(sys)` from
`Systems/RTS/build_rts.jl`, then
`build_model_with_flow_canceling_terms(sys; device_models = RTS_FC_DEVICE_MODELS)`.
