const DEFAULT_UC_MODELS = Dict(
    Line => StaticBranchBounds,
    TapTransformer => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicUnitCommitment,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    HydroDispatch => HydroDispatchRunOfRiver,
)

const DEFAULT_ED_MODELS = Dict(
    Line => StaticBranchBounds,
    TapTransformer => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicDispatch,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    HydroDispatch => HydroDispatchRunOfRiver,
)

const DEFAULT_MILP_OPTIMIZER = optimizer_with_attributes(Xpress.Optimizer)

const DEFAULT_NLP_OPTIMIZER = optimizer_with_attributes(Ipopt.Optimizer)

# PowerSimulations alias for convenience
const PSI = PowerSimulations

"""
Variable type representing the total approximated line losses in the system.
Used in optimization models to track aggregate transmission losses.
"""
struct LineLossTotalApproximation <: PSI.VariableType end

"""
Constraint type for approximating line loss relationships in the power system.
Links individual bus injections to total system losses via loss factors.
"""
struct LineLossConstraintApproximation <: PSI.ConstraintType end

"""
    make_ptdf_model_without_losses(
        sys::PSY.System
    ) -> DecisionModel

Create a PTDF-based decision model without building it.

# Arguments
- `sys::PSY.System`: PowerSystems.jl system containing network and device data
- `device_models`: Dictionary mapping device types to formulation models (default: DEFAULT_UC_MODELS)
- `optimizer`: JuMP-compatible optimizer (default: DEFAULT_MILP_OPTIMIZER)
- `ptdf`: Pre-computed PTDF matrix (default: nothing, computed if not provided)
- `name`: Name identifier for the model (default: "UC")

# Returns
- `DecisionModel`: Unbuilt decision model (call `build!` separately)

# Details
Creates the model structure with:
- PTDF-based network formulation
- Slack variables for constraint relaxation diagnostics
- Dual variable computation for nodal balance
- AC power flow evaluation with loss factor calculation

This function is used internally by `build_ptdf_model_without_losses`, which
additionally calls `build!` to construct the JuMP model.
"""
function make_ptdf_model_without_losses(
    sys::PSY.System;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
    ptdf = nothing,
    name = "UC",
    ignore_pf = false,
)

    # Create problem template with PTDF network model
    # - use_slacks: Add slack variables for infeasibility diagnosis
    # - duals: Compute dual variables for nodal balance constraints
    # - calculate_loss_factors: Enable loss factor computation in AC power flow
    if isnothing(ptdf)
        ptdf_used = PTDF(sys)
    else
        ptdf_used = ptdf
    end
    if !ignore_pf
        template_uc = ProblemTemplate(
            NetworkModel(
                PTDFPowerModel;
                PTDF_matrix = ptdf_used,
                use_slacks = true,
                duals = [CopperPlateBalanceConstraint],
                power_flow_evaluation = PowerFlows.ACPowerFlow(; calculate_loss_factors = true),
            ),
        )
    else
        template_uc = ProblemTemplate(
            NetworkModel(
                PTDFPowerModel;
                PTDF_matrix = ptdf_used,
                use_slacks = true,
                duals = [CopperPlateBalanceConstraint],
            ),
        )
    end

    # Set device models for all system components
    for (device_type, model) in device_models
        set_device_model!(template_uc, device_type, model)
    end

    # Create the decision model with hourly resolution
    decision_model = DecisionModel(
        template_uc,
        sys;
        optimizer = optimizer,
        name = name,
        store_variable_names = true,  # Store names for debugging
    )

    return decision_model
end

function make_acopf_model(
    sys::PSY.System;
    device_models = DEFAULT_ED_MODELS,
    optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf = nothing,
    name = "ED",
    ignore_pf = false,
)

    # Create problem template with PTDF network model
    # - use_slacks: Add slack variables for infeasibility diagnosis
    # - duals: Compute dual variables for nodal balance constraints
    # - calculate_loss_factors: Enable loss factor computation in AC power flow
    if isnothing(ptdf)
        ptdf_used = PTDF(sys)
    else
        ptdf_used = ptdf
    end
    if !ignore_pf
        template_uc = ProblemTemplate(
            NetworkModel(
                ACPPowerModel;
                use_slacks = true,
                power_flow_evaluation = PowerFlows.ACPowerFlow(; calculate_loss_factors = true),
            ),
        )
    else
        template_uc = ProblemTemplate(
            NetworkModel(
                ACPPowerModel;
                use_slacks = true,
            ),
        )
    end

    # Set device models for all system components
    for (device_type, model) in device_models
        set_device_model!(template_uc, device_type, model)
    end

    # Create the decision model with hourly resolution
    decision_model = DecisionModel(
        template_uc,
        sys;
        optimizer = optimizer,
        name = name,
        store_variable_names = true,  # Store names for debugging
    )

    return decision_model
end

"""
    build_ptdf_model_without_losses(
        sys::PSY.System,
        device_models = DEFAULT_UC_MODELS,
        optimizer = DEFAULT_MILP_OPTIMIZER
    ) -> DecisionModel

Build a PTDF (Power Transfer Distribution Factor) unit commitment model without losses.

# Arguments
- `sys::PSY.System`: PowerSystems.jl system containing network and device data
- `device_models`: Dictionary mapping device types to their formulation models (default: DEFAULT_UC_MODELS)
- `optimizer`: JuMP-compatible optimizer for solving the model (default: DEFAULT_MILP_OPTIMIZER)

# Returns
- `DecisionModel`: Built decision model ready for solving

# Details
Creates a NetworkModel with:
- PTDFPowerModel formulation
- Slack variables enabled for constraint relaxation
- Duals computed for CopperPlateBalanceConstraint
- AC power flow evaluation with loss factor calculation
"""
function build_ptdf_model_without_losses(
    sys::PSY.System;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
    ptdf = nothing,
    ignore_pf = false,
)
    decision_model = make_ptdf_model_without_losses(
        sys;
        device_models,
        optimizer,
        ptdf,
        ignore_pf,
    )
    build!(decision_model; output_dir = mktempdir())
    return decision_model
end

"""
    add_current_loss_variables!(model::PSI.DecisionModel) -> JuMPVariableArray

Add continuous decision variables for total line losses to the optimization model.

# Arguments
- `model::PSI.DecisionModel`: The decision model to modify

# Returns
- JuMP variable array indexed by reference bus and time step

# Details
Creates one loss variable per reference bus and time period to represent
the total system transmission losses at each time step.
"""
function add_current_loss_variables!(model::PSI.DecisionModel)
    # Access the optimization container
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)

    # Get reference buses from copper plate balance constraint
    con_bal = PSI.get_constraint(container, CopperPlateBalanceConstraint(), PSY.System)
    ref_buses = axes(con_bal, 1)

    # Create variable container for loss variables
    variable = PSI.add_variable_container!(
        container,
        LineLossTotalApproximation(),
        PSY.System,
        ref_buses,
        time_steps,
    )

    # Create JuMP variables for each reference bus and time step
    for ref_bus in ref_buses, t in time_steps
        variable[ref_bus, t] = JuMP.@variable(
            PSI.get_jump_model(container),
            base_name = "LineLossTotalApproximation_{$ref_bus}_{$t}"
        )
    end
    return variable
end

"""
    add_current_loss_to_copperplate_balance!(
        model::PSI.DecisionModel, 
        loss_variable, 
        total_loss_est
    )

Modify copper plate balance constraints to account for transmission losses.

# Arguments
- `model::PSI.DecisionModel`: The decision model to modify
- `loss_variable`: JuMP variable array representing total losses
- `total_loss_est`: Array of estimated total losses per time period (in MW)

# Details
Updates the nodal balance equation from (G - D = 0) to (G - D + loss_variable = total_loss_est)
where:
- G = total generation
- D = total demand
- loss_variable = optimization variable for losses
- total_loss_est = estimated losses from previous iteration

This formulation allows the optimizer to adjust generation to compensate for losses.
"""
function add_current_loss_to_copperplate_balance!(
    model::PSI.DecisionModel,
    loss_variable,
    total_loss_est,
)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)
    base_power = PSI.get_base_power(container)  # System base power for per-unit conversion

    # Get copper plate balance constraint (originally: G - D = 0)
    con_bal = PSI.get_constraint(container, CopperPlateBalanceConstraint(), PSY.System)
    ref_buses = axes(con_bal, 1)
    ref_bus = only(ref_buses)  # Ensure single reference bus exists

    for t in time_steps
        constraint = con_bal[ref_bus, t]

        # Add loss variable to left-hand side: (G - D) + loss_variable
        set_normalized_coefficient(constraint, loss_variable[ref_bus, t], 1)

        # Update right-hand side with loss estimate in per-unit
        # Final form: (G - D) + loss_variable = total_loss_est
        rhs = normalized_rhs(constraint)
        set_normalized_rhs(constraint, rhs + total_loss_est[t] / base_power)
    end
end

"""
    add_quadratic_current_loss_to_copperplate_balance!(
        model::PSI.DecisionModel, 
        loss_variable
    )

Modify copper plate balance constraints for quadratic loss formulation.

# Arguments
- `model::PSI.DecisionModel`: The decision model to modify
- `loss_variable`: JuMP variable array representing total losses

# Details
Updates the nodal balance equation from (G - D = 0) to (G - D + loss_variable = 0)
where:
- G = total generation
- D = total demand  
- loss_variable = optimization variable constrained by quadratic loss expression

**Difference from linear version:**
Unlike `add_current_loss_to_copperplate_balance!`, this function does NOT set
a fixed RHS value. Instead, the loss_variable will be constrained by a separate
quadratic constraint, allowing more accurate loss representation in NLP models.

The RHS remains zero because the quadratic loss constraint directly expresses
the loss_variable as a quadratic function of bus injections.
"""
function add_quadratic_current_loss_to_copperplate_balance!(
    model::PSI.DecisionModel,
    loss_variable,
)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)

    # Get copper plate balance constraint (originally: G - D = 0)
    con_bal = PSI.get_constraint(container, CopperPlateBalanceConstraint(), PSY.System)
    ref_buses = axes(con_bal, 1)
    ref_bus = only(ref_buses)  # Ensure single reference bus exists

    for t in time_steps
        constraint = con_bal[ref_bus, t]

        # Add loss variable to left-hand side: (G - D) + loss_variable = 0
        # The loss_variable value is determined by a separate quadratic constraint
        set_normalized_coefficient(constraint, loss_variable[ref_bus, t], 1)
    end
end

"""
    add_current_loss_constraint_approximation!(
        model, 
        loss_factors, 
        injection_old
    )

Add linearized loss approximation constraints to the optimization model.

# Arguments
- `model`: The decision model to modify
- `loss_factors`: Matrix of loss sensitivity factors (∂Loss/∂Injection) for each bus and time
- `injection_old`: Previous iteration's net injection values at each bus and time

# Details
Creates constraints linking total system losses to bus injections:

    loss_variable[t] = Σᵢ (injection_old[i,t] - injection[i,t]) * loss_factors[i,t]

This is a first-order Taylor expansion of losses around the previous operating point.
The constraint updates as the optimization adjusts bus injections from the previous values.
"""
function add_current_loss_constraint_approximation!(model, loss_factors, injection_old)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)

    # Retrieve the total loss variable
    loss_variable = PSI.get_variable(container, LineLossTotalApproximation(), PSY.System)
    ref_buses = axes(loss_variable, 1)

    # Create constraint container for loss approximation
    constraint = PSI.add_constraints_container!(
        container,
        LineLossConstraintApproximation(),
        PSY.System,
        ref_buses,
        time_steps,
    )

    # Get current iteration's injection expression
    injection = container.expressions[InfrastructureSystems.Optimization.ExpressionKey{
        ActivePowerBalance,
        ACBus,
    }(
        "",
    )]
    bus_ax = axes(injection, 1)
    num_buses = length(bus_ax)

    # Create linearized loss constraint for each time period
    # Loss ≈ Σ (change in injection) × (loss sensitivity at that bus)
    for ref_bus in ref_buses, t in time_steps
        constraint[ref_bus, t] = JuMP.@constraint(
            PSI.get_jump_model(container),
            loss_variable[ref_bus, t] == sum(
                (injection_old[i, t] - injection[bus_ax[i], t]) * loss_factors[i, t] for
                i in 1:num_buses
            )
        )
    end
end

"""
    add_current_loss_constraint_quadratic_approximation!(
        model, 
        sys, 
        ptdf, 
        res_old
    )

Add quadratic loss approximation constraints to the optimization model.

# Arguments
- `model`: The decision model to modify
- `sys`: PowerSystems.jl System object containing network data
- `ptdf`: PTDF matrix for the system
- `res_old`: Previous iteration's optimization results

# Details
Creates quadratic constraints linking total system losses to bus injections:

    loss_variable[t] = -Σₖ Rₖ × (Σⱼ (Vₖ,ₜ/Vⱼ,ₜ) × PTDFₖⱼ × Injectionⱼ,ₜ)²

where:
- Rₖ = resistance of branch k
- Vₖ,ₜ = voltage magnitude at branch k endpoint (from previous solution)
- Vⱼ,ₜ = voltage magnitude at bus j (from previous solution)
- PTDFₖⱼ = power transfer distribution factor
- Injectionⱼ,ₜ = net injection at bus j (decision variable)

**Formulation:**
This represents P = I²R losses where:
- Branch currents are approximated using PTDF and injections
- Voltage magnitudes are fixed from the previous AC power flow solution
- Results in a quadratic expression suitable for NLP solvers

**Advantages over linear:**
- More physically accurate loss representation
- Captures nonlinear relationship between power flow and losses
- Better suited for economic dispatch with fixed commitments
"""
function add_current_loss_constraint_quadratic_approximation!(model, sys, ptdf, res_old)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)
    arcs_length = length(axes(ptdf, 2))

    # Retrieve the total loss variable
    loss_variable = PSI.get_variable(container, LineLossTotalApproximation(), PSY.System)
    ref_buses = axes(loss_variable, 1)

    # Create constraint container for loss approximation
    constraint = PSI.add_constraints_container!(
        container,
        LineLossConstraintApproximation(),
        PSY.System,
        ref_buses,
        time_steps,
    )

    # Extract branch resistance values from system data
    R, _ = get_RX_vector(sys, ptdf)

    # Get current iteration's injection expression (decision variables)
    injection = container.expressions[InfrastructureSystems.Optimization.ExpressionKey{
        ActivePowerBalance,
        ACBus,
    }(
        "",
    )]
    bus_ax = axes(injection, 1)
    bus_length = length(bus_ax)

    # Extract voltage magnitudes from previous AC power flow solution
    V_bus = get_power_flow_voltage_mag(res_old)         # Bus voltages
    V_line = get_power_flow_arc_voltage_mag(res_old, ptdf)  # Branch endpoint voltages

    # Create quadratic loss constraint for each time period
    # Loss = Σ_branches R × (flow)² where flow is computed via PTDF
    for ref_bus in ref_buses, t in time_steps
        constraint[ref_bus, t] = JuMP.@constraint(
            PSI.get_jump_model(container),
            # Negative sign because losses reduce available power
            loss_variable[ref_bus, t] ==
            -sum(
                R[k] *
                (sum(
                    V_line[k, t] / V_bus[j, t] * ptdf[k, j] * injection[bus_ax[j], t]
                    for j in 1:bus_length
                ))^2 for k in 1:arcs_length
            )
        )
    end
end

"""
    update_copperplate_loss_approximation!(
        model::PSI.DecisionModel, 
        loss_factors, 
        total_loss_est
    )

Update the optimization model with linearized loss approximation for the current iteration.

# Arguments
- `model::PSI.DecisionModel`: The decision model to update
- `loss_factors`: Matrix of loss sensitivity factors for each bus and time
- `total_loss_est`: Array of estimated total losses per time period
- `injection_old`: Previous iteration's net injection values at each bus and time

# Details
This is a convenience function that orchestrates the three-step process:
1. Add loss variables to the model
2. Update copper plate balance to include losses
3. Add linearized constraints relating injections to losses
"""
function update_copperplate_loss_approximation!(
    model::PSI.DecisionModel,
    loss_factors,
    total_loss_est,
    injection_old,
)
    # Step 1: Add decision variables for total losses
    loss_var = add_current_loss_variables!(model)

    # Step 2: Modify nodal balance to account for losses
    add_current_loss_to_copperplate_balance!(model, loss_var, total_loss_est)

    # Step 3: Add linearized loss approximation constraints
    add_current_loss_constraint_approximation!(model, loss_factors, injection_old)
end

"""
    update_copperplate_quadratic_loss_approximation!(
        model::PSI.DecisionModel, 
        sys, 
        ptdf, 
        res_old
    )

Update the optimization model with quadratic loss approximation.

# Arguments
- `model::PSI.DecisionModel`: The decision model to update
- `sys`: PowerSystems.jl System object
- `ptdf`: PTDF matrix for the system
- `res_old`: Previous iteration's optimization results

# Details
This is a convenience function that orchestrates the three-step process
for adding quadratic losses:

1. Add loss variables to the model
2. Update copper plate balance to include loss variables
3. Add quadratic constraints relating bus injections to losses via P = I²R

**Use Case:**
Typically used for economic dispatch (ED) problems where:
- More accurate loss representation is needed
- NLP solver is available (quadratic constraints)
- Unit commitments are fixed from UC stage

**Comparison with linear:**
`update_copperplate_loss_approximation!` uses first-order Taylor approximation,
suitable for MILP. This function uses full quadratic formulation, more accurate
but requires NLP solver.
"""
function update_copperplate_quadratic_loss_approximation!(
    model::PSI.DecisionModel,
    sys,
    ptdf,
    res_old,
)
    # Step 1: Add decision variables for total losses
    loss_var = add_current_loss_variables!(model)

    # Step 2: Modify nodal balance to account for losses (without fixed RHS)
    add_quadratic_current_loss_to_copperplate_balance!(model, loss_var)

    # Step 3: Add quadratic loss approximation constraints (P = I²R formulation)
    add_current_loss_constraint_quadratic_approximation!(model, sys, ptdf, res_old)
end

"""
    update_transmission_constraints_with_losses!(
        model,
        res_old,
        sys
    )

Update PTDF network flow constraints to account for transmission losses.

# Arguments
- `model`: The decision model to modify
- `res_old`: Previous iteration's optimization results
- `sys`: PowerSystems.jl System object

# Details
This function modifies the right-hand side of branch flow constraints to incorporate
the effects of transmission losses using the Fictitious Nodal Demand (FND) approach.

**Method:**
1. Compute FND by distributing branch losses to endpoint buses
2. For each branch constraint, update RHS using PTDF matrix:
   
   Flow_limit_new = Flow_limit_old + Σᵢ PTDF[k,i] × FND[i,t]

3. Apply corrections to both Line and TapTransformer constraints

**Background:**
In a lossless PTDF model, branch flows are linear functions of injections.
To account for losses, we add fictitious demands at buses representing the
energy dissipated in nearby branches. The PTDF matrix then propagates these
fictitious demands through the network, adjusting branch flow limits accordingly.

This ensures that the transmission constraints reflect the additional loading
caused by losses, preventing constraint violations in the actual AC power flow.
"""
function update_transmission_constraints_with_losses!(model, res_old, sys, ptdf)
    container = model.internal.container
    base_power = PSI.get_base_power(container)  # System base power for per-unit conversion

    # Compute fictitious nodal demands representing distributed losses
    FND = get_fictitious_nodal_demand_by_loss(res_old, sys)

    # Get PTDF matrix and system topology information
    bus_axes = axes(ptdf, 1)
    num_bus = length(bus_axes)
    arc_axes = axes(ptdf, 2)

    # Get network flow constraints for lines and transformers
    con_line = PSI.get_constraint(container, NetworkFlowConstraint(), Line)
    con_tap = PSI.get_constraint(container, NetworkFlowConstraint(), TapTransformer)

    # Update line flow constraints
    for k in keys(con_line)
        rhs = normalized_rhs(con_line[k])
        line_name = k[1]
        time = k[2]

        # Find corresponding arc in PTDF matrix
        arc_ax = get_arc_axis_from_branch_name(sys, line_name)
        ptdf_ix = findfirst(x -> x == arc_ax, arc_axes)

        # Update RHS: add contribution from fictitious demands via PTDF
        set_normalized_rhs(
            con_line[k],
            rhs + sum(ptdf[ptdf_ix, i] * FND[i, time] / base_power for i in 1:num_bus),
        )
    end

    # Update tap transformer flow constraints
    for k in keys(con_tap)
        rhs = normalized_rhs(con_tap[k])
        tap_name = k[1]
        time = k[2]

        # Find corresponding arc in PTDF matrix
        arc_ax = get_arc_axis_from_branch_name(sys, tap_name)
        ptdf_ix = findfirst(x -> x == arc_ax, arc_axes)

        # Update RHS: add contribution from fictitious demands via PTDF
        set_normalized_rhs(
            con_tap[k],
            rhs + sum(ptdf[ptdf_ix, i] * FND[i, time] / base_power for i in 1:num_bus),
        )
    end
    return
end

"""
    build_ptdf_model_with_linear_losses(
        sys::PSY.System,
        res_old,
        injection_old
    ) -> DecisionModel

Build a PTDF unit commitment model with linearized transmission loss approximation.

# Arguments
- `sys::PSY.System`: PowerSystems.jl System object
- `res_old`: Previous iteration's optimization results (or lossless solution for first iteration)
- `injection_old`: Previous iteration's bus injection values

# Returns
- `DecisionModel`: Built decision model with loss approximation, ready for solving

# Details
This function creates a comprehensive loss-aware UC model by:

**Step 1:** Build baseline PTDF model (same structure as lossless model)

**Step 2:** Extract loss parameters from previous solution:
- `loss_factors`: ∂Loss/∂Injection sensitivity at each bus (from power flow)
- `total_loss_est`: Total system losses in MW (from AC power flow analysis)

**Step 3:** Update copper plate balance constraints:
- Add loss variables and constraints
- Linearize losses around previous operating point: 
  Loss ≈ Loss₀ + Σᵢ (∂Loss/∂Pᵢ) × (Pᵢ - Pᵢ,₀)

**Step 4:** Update transmission constraints:
- Adjust branch flow limits using Fictitious Nodal Demand (FND)
- Ensures consistency between nodal balance and branch flows
"""
function build_ptdf_model_with_linear_losses(
    sys::PSY.System,
    res_old,
    injection_old;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
)
    # Step 1: Build baseline PTDF model structure (without loss terms initially)
    model = build_ptdf_model_without_losses(sys; device_models, optimizer)

    # Step 2: Extract loss parameters from previous solution
    loss_factors = get_bus_loss_factors(res_old)      # ∂Loss/∂P at each bus
    total_loss_est = get_total_AC_loss(res_old)       # Total AC losses in MW

    # Step 3: Add linearized loss approximation to copper plate balance
    # This modifies the copperplate balance to account for generation needed to cover losses
    update_copperplate_loss_approximation!(
        model,
        loss_factors,
        total_loss_est,
        injection_old,
    )

    # Step 4: Update transmission constraints with fictitious nodal demands
    # This ensures branch flows reflect the additional loading from losses
    ptdf = PTDF(sys)
    update_transmission_constraints_with_losses!(model, res_old, sys, ptdf)

    return model
end

"""
    build_ptdf_model_with_quadratic_losses(
        sys::PSY.System,
        res_old
    ) -> DecisionModel

Build a PTDF economic dispatch model with quadratic transmission loss formulation.

# Arguments
- `sys::PSY.System`: PowerSystems.jl System object
- `res_old`: Previous iteration's optimization results (or UC solution for ED stage)
- `device_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `optimizer`: NLP optimizer for quadratic constraints (default: DEFAULT_NLP_OPTIMIZER)

# Returns
- `DecisionModel`: Built decision model with quadratic losses, ready for solving

# Details
This function creates an economic dispatch model with accurate quadratic loss representation:

**Step 1:** Build baseline PTDF model structure
- Uses ED device models (ThermalBasicDispatch, no commitment variables)
- Creates PTDF matrix for network representation

**Step 2:** Add quadratic loss approximation
- Losses modeled as: Loss = Σₖ Rₖ × (Iₖ)²
- Branch currents Iₖ computed via PTDF and bus injections
- Voltage magnitudes fixed from previous AC power flow
- Results in quadratic constraints requiring NLP solver

**Step 3:** Update transmission constraints
- Adjust branch flow limits using Fictitious Nodal Demand (FND)
- Ensures consistency between nodal balance and branch flows

**Typical Usage:**
Used as the Economic Dispatch stage in UC-ED simulations where:
- Unit commitments are fixed from UC solution (passed as res_old)
- More accurate loss modeling is desired for dispatch decisions
- NLP solver (Ipopt) can handle quadratic constraints

**Advantages over linear losses:**
- Physically accurate P = I²R relationship
- Better represents nonlinear loss behavior
- Improved dispatch optimality when losses are significant

# Example
```julia
# After solving UC with linear losses
uc_results = solve_uc_model(sys_uc)

# Build ED with quadratic losses using UC solution
ed_model = build_ptdf_model_with_quadratic_losses(sys_ed, uc_results)
solve!(ed_model)
```
"""
function build_ptdf_model_with_quadratic_losses(
    sys::PSY.System,
    res_old;
    device_models = DEFAULT_ED_MODELS,
    optimizer = DEFAULT_NLP_OPTIMIZER,
)
    # Compute PTDF matrix (needed for quadratic loss formulation)
    ptdf = PTDF(sys)

    # Step 1: Build baseline ED model structure (without loss terms initially)
    model = build_ptdf_model_without_losses(sys; device_models, optimizer, ptdf)

    # Step 2: Add quadratic loss formulation to copper plate balance
    # This creates Loss = Σ R × (PTDF × Injection)² constraints
    update_copperplate_quadratic_loss_approximation!(model, sys, ptdf, res_old)

    # Step 3: Update transmission constraints with fictitious nodal demands
    # Ensures branch flows reflect the additional loading from losses
    # update_transmission_constraints_with_losses!(model, res_old, sys, ptdf)

    return model
end
