"""
    run_uc_ed_lossless_simulation(
        sys_uc::PSY.System,
        sys_ed::PSY.System
    ) -> Tuple{Simulation, SimulationResults, Array, Array}

Run a multi-stage UC-ED simulation without transmission losses.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `uc_models`: Device models for UC stage (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED stage (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC problem (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED problem (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if not provided)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if not provided)

# Returns
- `sim`: Executed Simulation object
- `sim_res`: SimulationResults containing all stage results
- `injection_old_uc`: Array of net power injections at each bus for UC stage
- `injection_old_ed`: Array of net power injections at each bus for ED stage

# Details
This function creates and executes a sequential UC-ED simulation without loss modeling:
1. Builds a two-stage simulation with UC followed by ED
2. Uses PTDF formulation for both stages (lossless DC power flow)
3. Executes the simulation over the entire horizon
4. Extracts bus injection values from both stages for use in loss calculations

The injection values represent the starting point for iterative loss approximation
methods and can be used to initialize subsequent simulation runs with loss modeling.

# Example
```julia
sim, results, inj_uc, inj_ed = run_uc_ed_lossless_simulation(sys_uc, sys_ed)
```
"""
function run_uc_ed_lossless_simulation(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build the two-stage simulation (UC → ED) without loss constraints
    sim = build_uc_ed_simulation_with_no_losses(
        sys_uc,
        sys_ed;
        uc_models,
        ed_models,
        uc_optimizer,
        ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract internal containers for direct access to optimization models
    uc = sim.models.decision_models[1].internal.container
    ed = sim.models.decision_models[2].internal.container

    # Execute the full simulation sequence
    execute!(sim)

    # PSI clears JuMP solution values after execute!; calling optimize! directly on each
    # JuMP model restores them so that JuMP.value.(expr) returns correct values.
    optimize!(uc.JuMPmodel)
    optimize!(ed.JuMPmodel)

    # Extract bus injection values for both stages (Generation - Demand at each bus)
    injection_old_uc = deepcopy(
        JuMP.value.(
            uc.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
                ActivePowerBalance,
                ACBus,
            }(
                "",
            )]
        ).data,
    )
    injection_old_ed = deepcopy(
        JuMP.value.(
            ed.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
                ActivePowerBalance,
                ACBus,
            }(
                "",
            )]
        ).data,
    )

    # Package results for return
    sim_res = SimulationResults(sim)

    return sim, sim_res, injection_old_uc, injection_old_ed
end

"""
    run_uc_ed_acopf_simulation(
        sys_uc::PSY.System,
        sys_ed::PSY.System;
        ...
    ) -> Tuple{Simulation, SimulationResults, Array}

Run a cascaded UC+ED simulation where the ED stage uses an AC optimal power flow
formulation rather than a PTDF-based approximation.

The UC stage uses a lossless PTDF formulation (MILP). The ED stage solves an NLP
with full AC network constraints, providing accurate nodal voltages and branch flows.

# Arguments
- `sys_uc::PSY.System`: System for the Unit Commitment stage.
- `sys_ed::PSY.System`: System for the Economic Dispatch (ACOPF) stage.
- `uc_models`: Device models for UC (default: `DEFAULT_UC_MODELS`).
- `ed_models`: Device models for ED (default: `DEFAULT_ED_MODELS`).
- `uc_optimizer`: Optimizer for UC (default: `DEFAULT_MILP_OPTIMIZER`).
- `ed_optimizer`: Optimizer for ED NLP (default: `DEFAULT_NLP_OPTIMIZER`).
- `ptdf_uc`: Pre-computed PTDF matrix for UC (computed from `sys_uc` if `nothing`).
- `ptdf_ed`: Pre-computed PTDF matrix for ED (computed from `sys_ed` if `nothing`).

# Returns
- `(sim, sim_res, injection_old_uc)`: executed `Simulation`, `SimulationResults`,
  and UC-stage bus injection array (for use as a warm-start in iterative methods).
"""
function run_uc_ed_acopf_simulation(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build the two-stage simulation (UC → ED)
    sim = build_uc_ed_simulation_with_acopf(
        sys_uc,
        sys_ed;
        uc_models,
        ed_models,
        uc_optimizer,
        ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract internal containers for direct access to optimization models
    uc = sim.models.decision_models[1].internal.container
    ed = sim.models.decision_models[2].internal.container

    # Execute the full simulation sequence
    execute!(sim)

    # PSI clears JuMP solution values after execute!; calling optimize! directly on the
    # JuMP model restores them so that JuMP.value.(expr) returns correct values.
    optimize!(uc.JuMPmodel)

    # Extract bus injection values for both stages (Generation - Demand at each bus)
    injection_old_uc = deepcopy(
        JuMP.value.(
            uc.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
                ActivePowerBalance,
                ACBus,
            }(
                "",
            )]
        ).data,
    )

    # Package results for return
    sim_res = SimulationResults(sim)

    return sim, sim_res, injection_old_uc
end

"""
    run_uc_ed_quadratic_loss_simulation(
        sys_uc::PSY.System,
        sys_ed::PSY.System,
        res_old_uc,
        res_old_ed,
        injection_old_uc,
        injection_old_ed
    ) -> Tuple{Simulation, SimulationResults, Array, Array}

Run a multi-stage UC-ED simulation with hybrid loss modeling.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `res_old_uc`: Previous iteration's UC optimization results
- `res_old_ed`: Previous iteration's ED optimization results
- `injection_old_uc`: Previous iteration's UC bus injection values
- `injection_old_ed`: Previous iteration's ED bus injection values
- `uc_models`: Device models for UC stage (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED stage (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC problem (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED problem (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if not provided)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if not provided)

# Returns
- `sim`: Executed Simulation object with loss modeling
- `sim_res`: SimulationResults containing all stage results
- `injection_new_uc`: Updated array of net power injections from UC stage
- `injection_new_ed`: Updated array of net power injections from ED stage

# Details
This function implements a hybrid loss approximation strategy:
- **UC Stage**: Uses linearized transmission losses (MILP-compatible)
- **ED Stage**: Uses quadratic transmission losses (more accurate, requires NLP solver)

The loss terms are linearized/approximated around the operating point defined by
the previous iteration's injection values. This allows the optimization to account
for losses while maintaining tractability in the UC stage.

**Typical Usage Pattern:**
1. First run `run_uc_ed_lossless_simulation` to get initial operating point
2. Use those results as input to this function
3. Can be called iteratively until convergence

# Example
```julia
# First get lossless solution
sim0, res0, inj_uc0, inj_ed0 = run_uc_ed_lossless_simulation(sys_uc, sys_ed)
res_uc0 = get_decision_problem_results(res0, "UC")
res_ed0 = get_decision_problem_results(res0, "ED")

# Then run with losses
sim, res, inj_uc, inj_ed = run_uc_ed_quadratic_loss_simulation(
    sys_uc, sys_ed, res_uc0, res_ed0, inj_uc0, inj_ed0
)
```
"""
function run_uc_ed_quadratic_loss_simulation(
    sys_uc::PSY.System,
    sys_ed::PSY.System,
    res_old_uc,
    res_old_ed,
    injection_old_uc,
    injection_old_ed;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build simulation with hybrid loss modeling:
    # - UC uses linear losses (MILP-compatible)
    # - ED uses quadratic losses (more accurate)
    sim = build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses(
        sys_uc,
        sys_ed,
        res_old_uc,
        res_old_ed,
        injection_old_uc,
        injection_old_ed;
        uc_models = uc_models,
        ed_models = ed_models,
        uc_optimizer = uc_optimizer,
        ed_optimizer = ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract internal containers for direct access to optimization models
    uc = sim.models.decision_models[1].internal.container
    ed = sim.models.decision_models[2].internal.container

    # Execute the full simulation sequence with loss constraints
    execute!(sim)

    # PSI clears JuMP solution values after execute!; calling optimize! directly on each
    # JuMP model restores them so that JuMP.value.(expr) returns correct values.
    optimize!(uc.JuMPmodel)
    optimize!(ed.JuMPmodel)

    # Extract updated bus injection values for convergence checking
    injection_new_uc = deepcopy(
        JuMP.value.(
            uc.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
                ActivePowerBalance,
                ACBus,
            }(
                "",
            )]
        ).data,
    )
    injection_new_ed = deepcopy(
        JuMP.value.(
            ed.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
                ActivePowerBalance,
                ACBus,
            }(
                "",
            )]
        ).data,
    )

    # Package results for return
    sim_res = SimulationResults(sim)

    return sim, sim_res, injection_new_uc, injection_new_ed
end

"""
    run_iterative_uc_ed_quadratic_loss_simulation(
        sys_uc::PSY.System,
        sys_ed::PSY.System
    ) -> SimulationResults

Solve UC-ED simulation with transmission losses using successive approximations.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `uc_models`: Device models for UC stage (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED stage (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC problem (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED problem (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if not provided)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if not provided)
- `max_iter`: Maximum number of iterations (default: 10)
- `error_tol`: Convergence tolerance based on ED loss changes (default: 1e-3)

# Returns
- `SimulationResults`: Final converged simulation results

# Details
This function implements a successive approximation method for UC-ED with losses:

**Algorithm:**
1. Solve initial lossless UC-ED simulation (iteration 0)
2. For each iteration:
   a. Compute loss approximation around current operating point
   b. Solve UC-ED with:
      - UC: Linear loss approximation (MILP-compatible)
      - ED: Quadratic loss model (more accurate)
   c. Check convergence based on ED total loss changes
   d. If converged, return results
   e. Otherwise, update operating point and repeat

**Convergence:**
- Converges when absolute change in total ED losses < error_tol
- Monitors and prints diagnostic information each iteration:
  * Objective function changes
  * Total losses (time-series and sum)
  * Iteration error
- Stops early if converged or after max_iter iterations

**Loss Modeling Strategy:**
- UC uses linear losses to maintain MILP tractability for commitment decisions
- ED uses quadratic losses for more accurate dispatch with fixed commitments
- Each iteration refines approximation around new solution

**Performance Notes:**
- ED stage is focus of convergence metric (more sensitive to losses)
- PTDF matrices can be pre-computed and passed to avoid recomputation
- Typical convergence: 3-8 iterations for well-conditioned systems

# Example
```julia
# Basic usage with defaults
results = run_iterative_uc_ed_quadratic_loss_simulation(sys_uc, sys_ed)

# Custom convergence criteria
results = run_iterative_uc_ed_quadratic_loss_simulation(
    sys_uc, sys_ed;
    max_iter = 15,
    error_tol = 1e-4
)

# With pre-computed PTDF matrices
ptdf_uc = PTDF(sys_uc)
ptdf_ed = PTDF(sys_ed)
results = run_iterative_uc_ed_quadratic_loss_simulation(
    sys_uc, sys_ed;
    ptdf_uc = ptdf_uc,
    ptdf_ed = ptdf_ed
)
```
"""
function run_iterative_uc_ed_quadratic_loss_simulation(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
    max_iter = 10,
    error_tol = 1e-3,
)
    # Use provided PTDF matrices or compute them if not provided
    # Pre-computing these can save time in iterative runs
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Step 1: Solve initial lossless model to get starting point
    # This provides the base operating point for loss linearization
    sim_old, sim_res_old, injection_old_uc, injection_old_ed =
        run_uc_ed_lossless_simulation(
            sys_uc,
            sys_ed;
            uc_models = uc_models,
            ed_models = ed_models,
            uc_optimizer = uc_optimizer,
            ed_optimizer = uc_optimizer,  # Note: Using uc_optimizer for consistency
            ptdf_uc = ptdf_uc_used,
            ptdf_ed = ptdf_ed_used,
        )

    # Extract individual stage results for convergence tracking
    res_old_uc = get_decision_problem_results(sim_res_old, "UC")
    res_old_ed = get_decision_problem_results(sim_res_old, "ED")

    # Step 2: Iterate until convergence or max iterations
    for i in 1:max_iter
        println("Starting Iteration $i")

        # Solve simulation with losses approximated around previous solution
        sim_new, sim_res_new, injection_new_uc, injection_new_ed =
            run_uc_ed_quadratic_loss_simulation(
                sys_uc,
                sys_ed,
                res_old_uc,
                res_old_ed,
                injection_old_uc,
                injection_old_ed;
                uc_models = UC_MODELS,
                ed_models = ED_MODELS,
                uc_optimizer = uc_optimizer,
                ed_optimizer = ed_optimizer,
                ptdf_uc = ptdf_uc_used,
                ptdf_ed = ptdf_ed_used,
            )

        # Extract stage results for comparison
        res_new_uc = get_decision_problem_results(sim_res_new, "UC")
        res_new_ed = get_decision_problem_results(sim_res_new, "ED")

        # Compute convergence metrics and print diagnostics
        # Focus on ED stage as it has more accurate loss representation
        obj_func_old = read_optimizer_stats(res_old_ed)[1, "objective_value"]
        obj_func_new = read_optimizer_stats(res_new_ed)[1, "objective_value"]
        total_losses_old = get_total_AC_loss(res_old_ed)  # Time-series of losses
        total_losses_new = get_total_AC_loss(res_new_ed)

        # Print iteration diagnostics
        println("Current ED objective function difference: $(obj_func_new - obj_func_old)")
        println("Total ED loss previous iteration: $(total_losses_old)")
        println("Total ED loss current iteration: $(total_losses_new)")
        println("Total ED Loss Sum Previous Iteration: $(sum(total_losses_old))")
        println("Total ED Loss Sum Current Iteration: $(sum(total_losses_new))")

        # Convergence criterion: absolute change in total losses
        error_iteration = abs(sum(total_losses_new) - sum(total_losses_old))

        # Check for convergence
        if error_iteration < error_tol
            println("Finished iterative run at iteration $i with error $(error_iteration)")
            return sim_res_new
        end

        # Update operating point for next iteration
        res_old_uc = res_new_uc
        res_old_ed = res_new_ed
        injection_old_uc = injection_new_uc
        injection_old_ed = injection_new_ed

        # Handle max iterations reached without convergence
        if i == max_iter
            println(
                "Reached maximum number of iterations ($max_iter) with error $(error_iteration)",
            )
            return sim_res_new
        end
    end
end

"""
    run_uc_linear_loss_ed_acopf_simulation(
        sys_uc, sys_ed, res_old_uc, res_old_ed, injection_old_uc; ...
    ) -> Tuple{Simulation, SimulationResults, Array}

Run one UC+ED iteration using linearised losses in UC and an AC optimal power flow
in the ED stage, with both linearisations anchored to a previous operating point.

# Arguments
- `sys_uc::PSY.System`: System for the Unit Commitment stage.
- `sys_ed::PSY.System`: System for the Economic Dispatch (ACOPF) stage.
- `res_old_uc`: Previous iteration's UC `OptimizationProblemResults` (used to build
  the linear loss approximation for UC).
- `res_old_ed`: Previous iteration's ED results (used to anchor the ED linearisation).
- `injection_old_uc`: Previous iteration's UC bus injection array (p.u., all time-steps).
- `uc_models`, `ed_models`, `uc_optimizer`, `ed_optimizer`: see `run_uc_ed_lossless_simulation`.
- `ptdf_uc`, `ptdf_ed`: Pre-computed PTDF matrices (computed if `nothing`).

# Returns
- `(sim, sim_res, injection_new_uc)`: executed `Simulation`, `SimulationResults`,
  and updated UC-stage bus injection array for the next iteration.
"""
function run_uc_linear_loss_ed_acopf_simulation(
    sys_uc::PSY.System,
    sys_ed::PSY.System,
    res_old_uc,
    res_old_ed,
    injection_old_uc;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build simulation with hybrid loss modeling:
    # - UC uses linear losses (MILP-compatible)
    # - ED uses quadratic losses (more accurate)
    sim = build_uc_ed_simulation_with_acopf_and_uc_linear_losses(
        sys_uc,
        sys_ed,
        res_old_uc,
        res_old_ed,
        injection_old_uc;
        uc_models = uc_models,
        ed_models = ed_models,
        uc_optimizer = uc_optimizer,
        ed_optimizer = ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract internal containers for direct access to optimization models
    uc = sim.models.decision_models[1].internal.container
    ed = sim.models.decision_models[2].internal.container

    # Execute the full simulation sequence with loss constraints
    execute!(sim)

    # PSI clears JuMP solution values after execute!; calling optimize! directly on the
    # JuMP model restores them so that JuMP.value.(expr) returns correct values.
    optimize!(uc.JuMPmodel)

    # Extract updated bus injection values for convergence checking
    injection_new_uc = deepcopy(
        JuMP.value.(
            uc.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
                ActivePowerBalance,
                ACBus,
            }(
                "",
            )]
        ).data,
    )

    # Package results for return
    sim_res = SimulationResults(sim)

    return sim, sim_res, injection_new_uc
end

"""
    run_iterative_uc_linear_ed_acopf_simulation(
        sys_uc::PSY.System,
        sys_ed::PSY.System; ...
    ) -> SimulationResults

Solve a cascaded UC+ED simulation iteratively, using linearised losses in UC and
a full AC optimal power flow in the ED stage, until convergence or `max_iter`.

The initial operating point comes from `run_uc_ed_acopf_simulation` (lossless UC +
ACOPF ED). Each subsequent iteration calls `run_uc_linear_loss_ed_acopf_simulation`
with the previous solution as the linearisation anchor.

Convergence is declared when the absolute change in total ED AC losses between
successive iterations falls below `error_tol` — i.e., when re-linearising the
losses around the new operating point no longer meaningfully shifts the dispatch.

# Arguments
- `sys_uc::PSY.System`: System for the Unit Commitment stage.
- `sys_ed::PSY.System`: System for the Economic Dispatch (ACOPF) stage.
- `uc_models`, `ed_models`, `uc_optimizer`, `ed_optimizer`: see `run_uc_ed_lossless_simulation`.
- `ptdf_uc`, `ptdf_ed`: Pre-computed PTDF matrices (computed if `nothing`).
- `max_iter`: Maximum number of iterations (default: `5`).
- `error_tol`: Convergence threshold on absolute total-loss change (default: `1e-1` MW).

# Returns
- `SimulationResults` from the final (converged or last) iteration.
"""
function run_iterative_uc_linear_ed_acopf_simulation(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
    max_iter = 5,
    error_tol = 1e-1,
)
    # Use provided PTDF matrices or compute them if not provided
    # Pre-computing these can save time in iterative runs
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Step 1: Solve initial lossless model to get starting point
    # This provides the base operating point for loss linearization
    sim_old, sim_res_old, injection_old_uc =
        run_uc_ed_acopf_simulation(
            sys_uc,
            sys_ed;
            uc_models = uc_models,
            ed_models = ed_models,
            uc_optimizer = uc_optimizer,
            ed_optimizer = ed_optimizer,
            ptdf_uc = ptdf_uc_used,
            ptdf_ed = ptdf_ed_used,
        )

    # Extract individual stage results for convergence tracking
    res_old_uc = get_decision_problem_results(sim_res_old, "UC")
    res_old_ed = get_decision_problem_results(sim_res_old, "ED")

    # Step 2: Iterate until convergence or max iterations
    for i in 1:max_iter
        println("Starting Iteration $i")

        # Solve simulation with losses approximated around previous solution
        sim_new, sim_res_new, injection_new_uc =
            run_uc_linear_loss_ed_acopf_simulation(
                sys_uc,
                sys_ed,
                res_old_uc,
                res_old_ed,
                injection_old_uc;
                uc_models = UC_MODELS,
                ed_models = ED_MODELS,
                uc_optimizer = uc_optimizer,
                ed_optimizer = ed_optimizer,
                ptdf_uc = ptdf_uc_used,
                ptdf_ed = ptdf_ed_used,
            )

        # Extract stage results for comparison
        res_new_uc = get_decision_problem_results(sim_res_new, "UC")
        res_new_ed = get_decision_problem_results(sim_res_new, "ED")

        # Compute convergence metrics and print diagnostics
        # Focus on ED stage as it has more accurate loss representation
        obj_func_old = read_optimizer_stats(res_old_uc)[1, "objective_value"]
        obj_func_new = read_optimizer_stats(res_new_uc)[1, "objective_value"]
        total_losses_old = get_total_AC_loss_CATS(res_old_ed)  # Time-series of losses
        total_losses_new = get_total_AC_loss_CATS(res_new_ed)

        # Print iteration diagnostics
        println("Current UC objective function difference: $(obj_func_new - obj_func_old)")
        println("Total ED loss previous iteration: $(total_losses_old)")
        println("Total ED loss current iteration: $(total_losses_new)")
        println("Total ED Loss Sum Previous Iteration: $(sum(total_losses_old))")
        println("Total ED Loss Sum Current Iteration: $(sum(total_losses_new))")

        # Convergence criterion: absolute change in total losses
        error_iteration = abs(sum(total_losses_new) - sum(total_losses_old))

        # Check for convergence
        if error_iteration < error_tol
            println("Finished iterative run at iteration $i with error $(error_iteration)")
            return sim_res_new
        end

        # Update operating point for next iteration
        res_old_uc = res_new_uc
        res_old_ed = res_new_ed
        injection_old_uc = injection_new_uc

        # Handle max iterations reached without convergence
        if i == max_iter
            println(
                "Reached maximum number of iterations ($max_iter) with error $(error_iteration)",
            )
            return sim_res_new
        end
    end
end