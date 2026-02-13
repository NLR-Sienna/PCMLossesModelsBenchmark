"""
    run_lossless_model(sys::PSY.System) -> Tuple{DecisionModel, OptimizationProblemResults, Array}

Build and solve a unit commitment model without transmission losses.

# Arguments
- `sys::PSY.System`: PowerSystems.jl System object containing network and device data

# Returns
- `model`: Solved DecisionModel object
- `res`: OptimizationProblemResults containing solution data
- `injection_vals`: Array of net active power injections at each bus

# Details
This function creates a baseline solution assuming lossless transmission:
1. Builds a PTDF-based unit commitment model without loss terms
2. Solves the optimization problem
3. Extracts the results and bus injection values

The injection values represent the net power balance (Generation - Demand) at
each bus and are used as the starting point for iterative loss calculations.

This is typically used as the first step in iterative loss approximation methods.
"""
function run_lossless_model(
    sys::PSY.System;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER
)
    # Build the optimization model without loss considerations
    model = build_ptdf_model_without_losses(sys; device_models = device_models, optimizer = optimizer)
    
    # Solve the optimization problem
    solve!(model)

    # Extract results
    res = OptimizationProblemResults(model)
    
    # Re-optimize to ensure fresh solution (workaround for result extraction)
    JuMP.optimize!(model.internal.container.JuMPmodel)
    
    # Extract net injection values at each bus for loss calculation
    injection_vals = deepcopy(JuMP.value.(model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")]).data)
    
    return model, res, injection_vals
end

"""
    run_linear_loss_model(
        sys::PSY.System,
        res_old,
        injection_old
    ) -> Tuple{DecisionModel, OptimizationProblemResults, Array}

Build and solve a unit commitment model with linearized transmission losses.

# Arguments
- `sys::PSY.System`: PowerSystems.jl System object
- `res_old`: Previous iteration's optimization results
- `injection_old`: Previous iteration's bus injection values

# Returns
- `model`: Solved DecisionModel with loss approximation
- `res`: OptimizationProblemResults containing solution data
- `injection_vals`: Updated array of net active power injections

# Details
This function creates an improved solution accounting for transmission losses:
1. Builds a PTDF model with linearized loss approximation based on previous solution
2. Solves the optimization problem with updated loss terms
3. Extracts results and injection values for potential next iteration
"""
function run_linear_loss_model(
    sys::PSY.System,
    res_old,
    injection_old;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER
)
    # Build model with linear loss approximation around previous operating point
    model = build_ptdf_model_with_linear_losses(sys, res_old, injection_old; device_models = device_models, optimizer = optimizer)
    
    # Solve the optimization problem
    solve!(model)

    # Extract results
    res = OptimizationProblemResults(model)
    
    # Re-optimize to ensure fresh solution (workaround for result extraction)
    JuMP.optimize!(model.internal.container.JuMPmodel)
    
    # Extract updated injection values for convergence checking and next iteration
    injection_vals = deepcopy(JuMP.value.(model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")]).data)
    
    return model, res, injection_vals
end

"""
    run_iterative_linear_loss_model(
        sys::PSY.System,
        max_iter = 15,
        error_tol = 1e-3
    ) -> OptimizationProblemResults

Solve unit commitment with transmission losses using successive linear approximations.

# Arguments
- `sys::PSY.System`: PowerSystems.jl System object
- `max_iter`: Maximum number of iterations (default: 15)
- `error_tol`: Convergence tolerance based on generator output changes (default: 1e-3)

# Returns
- `OptimizationProblemResults`: Final converged solution results

# Details
This function implements a successive linearization method for unit commitment with losses:

**Algorithm:**
1. Solve initial lossless model (iteration 0)
2. For each iteration:
   a. Linearize losses around current operating point
   b. Solve UC with linear loss approximation
   c. Check convergence based on generator output changes
   d. If converged, return results
   e. Otherwise, update operating point and repeat

**Convergence:**
- Converges when generator output changes between iterations < error_tol
- Error computed via `compute_iterative_error_based_on_generator_output`
- Stops early if converged or after max_iter iterations

**Benefits:**
- More accurate than single-pass linear approximation
- Accounts for interaction between dispatch and losses
- Each iteration refines the loss approximation around new solution

# Example
```julia
res = run_iterative_linear_loss_model(sys, 20, 1e-4)  # Custom tolerances
```
"""
function run_iterative_linear_loss_model(
    sys::PSY.System;
    max_iter = 15,
    error_tol = 1e-3,
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
)
    # Step 1: Solve initial lossless model to get starting point
    model_old, res_old, injection_old = run_lossless_model(sys; device_models = device_models, optimizer = optimizer);
    
    # Step 2: Iterate until convergence or max iterations
    for i = 1:max_iter
        println("Starting Iteration $i")
        
        # Solve model with losses linearized around previous solution
        model_new, res_new, injection_new = run_linear_loss_model(sys, res_old, injection_old; device_models = device_models, optimizer = optimizer);
        
        # Compute convergence metric (change in generator outputs)
        error_iteration = compute_iterative_error_based_on_generator_output(res_old, res_new)
        println("Current error: $(error_iteration)")
        obj_func_old = res_old.optimizer_stats[1, "objective_value"]
        obj_func_new = res_new.optimizer_stats[1, "objective_value"]
        total_losses_old = get_total_AC_loss(res_old)
        total_losses_new = get_total_AC_loss(res_new)
        println("Current objective function difference: $(obj_func_new - obj_func_old)")
        println("Total loss previous iteration: $(total_losses_old)")
        println("Total loss current iteration: $(total_losses_new)")
        println("Total Loss Sum Previous Iteration: $(sum(total_losses_old))")
        println("Total Loss Sum Current Iteration: $(sum(total_losses_new))")
        
        # Check for convergence
        if error_iteration < error_tol
            println("Finished iterative run at iteration $i with error $(error_iteration)")
            return res_new
        end
        
        # Update for next iteration
        res_old = res_new
        injection_old = injection_new
        
        # Handle max iterations reached
        if i == max_iter
            println("Reached maximum number of iterations ($max_iter) with error $(error_iteration)")
            return res_new
        end
    end
end
