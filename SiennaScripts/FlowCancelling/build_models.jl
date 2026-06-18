const DEFAULT_UC_MODELS = Dict(
    Line => StaticBranchBounds,
    #TapTransformer => StaticBranch,
    #Transformer2W => StaticBranch,
    PhaseShiftingTransformer => StaticBranch,
    ThermalStandard => ThermalDispatchNoMin,
    PowerLoad => StaticPowerLoad,
    #RenewableDispatch => RenewableFullDispatch,
    #TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    #HydroDispatch => HydroDispatchRunOfRiver,
)

# Xpress is the intended MILP solver. Guard the reference so the file still loads in
# environments where Xpress is unavailable (e.g. CI/testing); callers can then pass an
# `optimizer` keyword (e.g. HiGHS) to the build functions instead.
const DEFAULT_MILP_OPTIMIZER = try
    optimizer_with_attributes(Xpress.Optimizer)
catch
    nothing
end

import PowerNetworkMatrices

const PSI = PowerSimulations
const PSY = PowerSystems
const M_max = 10.0

# Default device-model dicts used to build the flow-cancelling templates.
# Branches must use `StaticBranch` so the `FlowRateConstraint` containers (which the
# flow-cancelling code mutates) are created.
const DEFAULT_FC_DEVICE_MODELS = Dict(
    ThermalStandard => ThermalDispatchNoMin,
    Line => StaticBranch,
    PhaseShiftingTransformer => StaticBranch,
    PowerLoad => StaticPowerLoad,
)

# RTS carries TapTransformer plus renewables/hydro/HVDC, so it needs a richer dict.
const RTS_FC_DEVICE_MODELS = Dict(
    ThermalStandard => ThermalDispatchNoMin,
    Line => StaticBranch,
    TapTransformer => StaticBranch,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalDispatch,
)

struct GenerationInvestmentVariable <: PSI.VariableType end
struct BranchInvestmentVariable <: PSI.VariableType end
struct BranchCancellingFlowVariable <: PSI.VariableType end
struct BigMConstraint <: PSI.ConstraintType end


"""
    make_base_ptdf_model(sys; device_models, optimizer, ptdf, name, ignore_pf) -> DecisionModel

Create a PTDF-based `DecisionModel` without building it.  Use `build_base_ptdf_model` to
also call `build!`.  When `ignore_pf = false` the template includes an AC power flow
evaluation with loss factor calculation.
"""
function make_base_ptdf_model(
    sys::PSY.System;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
    ptdf = nothing,
    name = "UC",
    ignore_pf = true,
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
                power_flow_evaluation = PowerFlows.ACPowerFlow(; calculate_loss_factors = true, calculate_voltage_stability_factors = true,
            )),
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

"""
    build_base_ptdf_model(sys; device_models, optimizer, ptdf, ignore_pf) -> DecisionModel

Build and return a PTDF-based `DecisionModel` ready for post-processing.  Internally calls
`make_base_ptdf_model` followed by `build!`.
"""
function build_base_ptdf_model(
    sys::PSY.System;
    device_models = DEFAULT_UC_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
    ptdf = nothing,
    ignore_pf = true,
)
    decision_model = make_base_ptdf_model(
        sys;
        device_models,
        optimizer,
        ptdf,
        ignore_pf,
    )
    build!(decision_model; output_dir = mktempdir())
    return decision_model
end

# ---------------------------------------------------------------------------
# Flow-cancelling (Caramanis et al., IEEE PES 2016) post-processing
# ---------------------------------------------------------------------------

"""
    _branch_arc_map(ptdf) -> name_to_arc_map

Populate (if needed) and return the network-reduction `name_to_arc_map`, a
`Dict{DataType, Dict{String, Tuple{Tuple{Int,Int}, String}}}` mapping each branch type to
a `reduced-name -> (arc, source)` dictionary.  The *reduced* names are exactly the names
used on the `FlowRateConstraint` / `PTDFBranchFlow` axes — including the `-double_circuit`
names produced when parallel branches are aggregated.  This is the bridge between
constraint-axis names and PTDF rows for systems with parallel lines.
"""
function _branch_arc_map(ptdf)
    PowerNetworkMatrices.populate_branch_maps_by_type!(ptdf.network_reduction_data)
    return ptdf.network_reduction_data.name_to_arc_map
end

"""
    _resolve_branch_arc(name_to_arc_map, sys, T, name) -> Tuple{Int,Int}

Resolve the `(from_bus_number, to_bus_number)` arc for a branch named `name` (of type `T`)
as it appears on a constraint/expression axis.  Handles aggregated `-double_circuit` names
via `name_to_arc_map`; falls back to a component lookup (stripping the suffix) if a name is
unexpectedly absent from the map.
"""
function _resolve_branch_arc(name_to_arc_map, sys, T, name::AbstractString)
    if haskey(name_to_arc_map, T) && haskey(name_to_arc_map[T], name)
        return first(name_to_arc_map[T][name])  # (arc, source) -> arc
    end
    @warn "Branch `$name` (type $T) not in name_to_arc_map; falling back to component lookup."
    return get_arc_axis_from_branch_name(sys, name)
end

"""
    _get_ptdf_shift_term(ptdf, monitored_arc, from_bus_num, to_bus_num) -> Float64

Return the shift factor for the branch whose arc is `monitored_arc` (a `(from,to)` bus-number
tuple, i.e. the PTDF row key) due to injecting at `from_bus_num` minus the factor due to
injecting at `to_bus_num` (eq. 3: Δ_{l,k} = H_{from_k,l} − H_{to_k,l}).
Returns 0.0 with a warning when the arc or a bus is not found in `ptdf`.

`monitored_arc` is resolved from a constraint-axis name via [`_resolve_branch_arc`], which
keeps the lookup correct for aggregated `-double_circuit` parallel branches.
"""
function _get_ptdf_shift_term(
    ptdf,
    monitored_arc::Tuple{Int,Int},
    from_bus_num::Int,
    to_bus_num::Int,
)
    bus_lookup = ptdf.lookup[1]   # bus number  -> col index
    br_lookup  = ptdf.lookup[2]   # branch arc  -> row index

    if !haskey(br_lookup, monitored_arc)
        @warn "Monitored branch arc $monitored_arc not found in PTDF; shift term set to 0."
        return 0.0
    end
    if !haskey(bus_lookup, from_bus_num) || !haskey(bus_lookup, to_bus_num)
        @warn "Bus $from_bus_num or $to_bus_num not found in PTDF; shift term set to 0."
        return 0.0
    end

    line_arc_ix = br_lookup[monitored_arc]
    from_bus_ix = bus_lookup[from_bus_num]
    to_bus_ix   = bus_lookup[to_bus_num]
    return ptdf.data[from_bus_ix, line_arc_ix] - ptdf.data[to_bus_ix, line_arc_ix]
end

"""
    add_branch_investment_variables!(decision_model, candidate_lines, T) -> variable container

Add one binary investment variable `z_k ∈ {0,1}` per entry in `candidate_lines` to the PSI
container and the underlying JuMP model, keyed by component type `T`.  Returns the variable
container indexed by candidate line name.
"""
function add_branch_investment_variables!(decision_model, candidate_lines, T)
    container  = decision_model.internal.container
    jump_model = PSI.get_jump_model(container)
    names      = get_name.(candidate_lines)

    variable = PSI.add_variable_container!(
        container,
        BranchInvestmentVariable(),
        T,
        names,
    )

    for name in names
        variable[name] = JuMP.@variable(
            jump_model,
            binary    = true,
            base_name = "BranchInvestment_$(T)_{$name}",
        )
    end
    return variable
end

"""
    add_branch_cancelling_flow_variables!(decision_model, candidate_lines, T) -> variable container

Add a continuous flow-cancelling variable `v_k[t]` per candidate line and time step to the
PSI container, keyed by component type `T`.  `v_k[t] ≈ z_k · FlowVar[k,t]` (enforced via
big-M constraints).  Returns the variable container indexed by `(candidate_line_name, t)`.
"""
function add_branch_cancelling_flow_variables!(decision_model, candidate_lines, T)
    container  = decision_model.internal.container
    jump_model = PSI.get_jump_model(container)
    names      = get_name.(candidate_lines)
    time_steps = PSI.get_time_steps(container)

    variable = PSI.add_variable_container!(
        container,
        BranchCancellingFlowVariable(),
        T,
        names,
        time_steps,
    )

    for name in names, t in time_steps
        variable[name, t] = JuMP.@variable(
            jump_model,
            base_name = "BranchCancellingFlow_$(T)_{$name, $t}",
        )
    end
    return variable
end

"""
    add_bigM_linking_constraints!(decision_model, candidate_lines, T, z_var, v_var)

Add big-M constraints linking the flow-cancelling variable `v_k[t]` to the binary investment
variable `z_k`, keyed by component type `T`:

    v_k[t] ≤  M · (1 − z_k)
    v_k[t] ≥ −M · (1 − z_k)

where M = `M_max` (global constant).  When `z_k = 1` (line built) `v_k` is forced to zero;
when `z_k = 0` (line not built) `v_k` is free within ±M.
"""
function add_bigM_linking_constraints!(decision_model, candidate_lines, T, z_var, v_var)
    container  = decision_model.internal.container
    jump_model = PSI.get_jump_model(container)
    time_steps = PSI.get_time_steps(container)

    constraint_ub =
        PSI.add_constraints_container!(
            container,
            BigMConstraint(),
            T,
            get_name.(candidate_lines),
            time_steps;
            meta = "ub"
        )
    constraint_lb =
        PSI.add_constraints_container!(
            container,
            BigMConstraint(),
            T,
            get_name.(candidate_lines),
            time_steps;
            meta = "lb"
        )

    for line in candidate_lines
        name = get_name(line)
        # v_k is bounded by the candidate's own rating (its flow when built), so a
        # rating-based big-M is the natural, non-binding choice. RTS flows exceed the old
        # fixed M_max = 10 p.u., which would have artificially constrained v_k.
        Big_M = max(get_rating(line), M_max)
        for t in time_steps
            constraint_ub[name, t] = JuMP.@constraint(jump_model, v_var[name, t] <= Big_M * (1 - z_var[name]))
            constraint_lb[name, t] = JuMP.@constraint(jump_model, v_var[name, t] >= -Big_M * (1 - z_var[name]))
        end
    end
    return
end

"""
    add_shift_terms_to_existing_line_constraints!(
        decision_model, T, candidate_lines, v_var, ptdf, sys, name_to_arc_map)

For each existing (non-candidate) branch of type `T` and each candidate line `k`, append
`Δ_{l,k} · v_k[t]` to the `FlowRateConstraint` upper and lower bounds (eq. 3).
`Δ_{l,k} = PTDF[from_k, l] − PTDF[to_k, l]`.

Iteration is over the **constraint axes** (the reduced branch names used by the model),
not over PSY components, so aggregated `-double_circuit` parallel branches are handled
correctly: each parallel pair has a single reduced constraint and is visited exactly once.
Candidate names are skipped here (their constraints are modified by
[`add_shift_terms_to_candidate_line_constraints!`]).  Returns immediately if `T` has no
`FlowRateConstraint` container in the model.
"""
function add_shift_terms_to_existing_line_constraints!(
    decision_model,
    T,
    candidate_lines,
    v_var,
    ptdf,
    sys,
    name_to_arc_map,
)
    container  = decision_model.internal.container
    time_steps = PSI.get_time_steps(container)

    ub_key = InfrastructureSystems.Optimization.ConstraintKey{FlowRateConstraint, T}("ub")
    lb_key = InfrastructureSystems.Optimization.ConstraintKey{FlowRateConstraint, T}("lb")
    (haskey(container.constraints, ub_key) && haskey(container.constraints, lb_key)) ||
        return nothing
    ub_con = container.constraints[ub_key]
    lb_con = container.constraints[lb_key]

    candidate_names = Set(get_name.(candidate_lines))
    cand_info = [
        (get_name(c), get_number(get_from(get_arc(c))), get_number(get_to(get_arc(c))))
        for c in candidate_lines
    ]

    for l_name in axes(ub_con, 1)
        l_name in candidate_names && continue
        monitored_arc = _resolve_branch_arc(name_to_arc_map, sys, T, l_name)
        for (k_name, from_num, to_num) in cand_info
            shift = _get_ptdf_shift_term(ptdf, monitored_arc, from_num, to_num)
            for t in time_steps
                set_normalized_coefficient(ub_con[l_name, t], v_var[k_name, t], shift)
                set_normalized_coefficient(lb_con[l_name, t], v_var[k_name, t], shift)
            end
        end
    end
    return nothing
end

"""
    add_shift_terms_to_candidate_line_constraints!(
        decision_model, sys, candidate_lines, T, z_var, v_var, ptdf)

Modify the `FlowRateConstraint` bounds for each candidate line to implement eq. (18):

    UB:  FlowVar[k,t] − rating_k · z_k + (Δ_{k,k} − 1) · v_k[t] + Σ_{j≠k} Δ_{k,j} · v_j[t] ≤ 0
    LB: −FlowVar[k,t] − rating_k · z_k − (Δ_{k,k} − 1) · v_k[t] − Σ_{j≠k} Δ_{k,j} · v_j[t] ≤ 0

The RHS is shifted from ±rating_k to 0.  `z_k = 0` forces `FlowVar[k,t] = 0`;
`z_k = 1` restores a ±rating_k feasible range via the `v_k` cancellation.
"""
function add_shift_terms_to_candidate_line_constraints!(
    decision_model,
    sys,
    candidate_lines,
    T,
    z_var,
    v_var,
    ptdf,
    name_to_arc_map,
)
    container  = decision_model.internal.container
    time_steps = PSI.get_time_steps(container)
    ub_con = PSI.get_constraint(container, FlowRateConstraint(), T, "ub")
    lb_con = PSI.get_constraint(container, FlowRateConstraint(), T, "lb")

    for line in candidate_lines
        k_name   = get_name(line)
        arc_k    = get_arc(line)
        k_from   = get_number(get_from(arc_k))
        k_to     = get_number(get_to(arc_k))
        rating_k = get_rating(line)

        # Candidate lines are never parallel (intermediate-bus trick), so the monitored
        # arc is just the candidate's own arc.
        monitored_arc_k = _resolve_branch_arc(name_to_arc_map, sys, T, k_name)

        # Self shift: Δ_{k,k} = H_{k, from_k} - H_{k, to_k}
        self_shift = _get_ptdf_shift_term(ptdf, monitored_arc_k, k_from, k_to)

        for t in time_steps
            ub = ub_con[k_name, t]
            lb = lb_con[k_name, t]

            # --- modify RHS: from ±rating_k to 0 ---
            ub_rhs = normalized_rhs(ub)
            set_normalized_rhs(ub, ub_rhs - rating_k)
            lb_rhs = normalized_rhs(lb)
            set_normalized_rhs(lb, lb_rhs + rating_k)

            # --- binary variable with coefficient -rating_k ---
            set_normalized_coefficient(ub, z_var[k_name], -rating_k)
            set_normalized_coefficient(lb, z_var[k_name], rating_k)

            # --- self flow-cancelling variable: (Δ_{k,k} - 1) for UB, -(Δ_{k,k}-1) for LB ---
            set_normalized_coefficient(ub, v_var[k_name, t],  (self_shift - 1.0))
            set_normalized_coefficient(lb, v_var[k_name, t], (self_shift - 1.0))

            # --- cross-shift terms from other candidate lines ---
            for other in candidate_lines
                kk_name = get_name(other)
                kk_name == k_name && continue

                arc_kk    = get_arc(other)
                kk_from   = get_number(get_from(arc_kk))
                kk_to     = get_number(get_to(arc_kk))
                cross_shift = _get_ptdf_shift_term(ptdf, monitored_arc_k, kk_from, kk_to)

                set_normalized_coefficient(ub, v_var[kk_name, t], cross_shift)
                set_normalized_coefficient(lb, v_var[kk_name, t], cross_shift)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Investment cost and generation investment constraints
# ---------------------------------------------------------------------------

struct GenerationInvestmentConstraint <: PSI.ConstraintType end

"""
    add_candidate_line_investment_costs!(decision_model, z_var)

For every candidate line (those whose name contains `"candidate"`), add
`project_cost * z_k` to the objective function, where `project_cost` is
read from `get_ext(line)["project_cost"]` (defaults to 0 if absent).

`z_var` is the binary-variable container returned by
`add_flow_canceling_terms!`.
"""
function add_candidate_line_investment_costs!(decision_model, z_var)
    sys        = PSI.get_system(decision_model)
    container  = decision_model.internal.container
    jump_model = PSI.get_jump_model(container)

    candidate_lines = filter(
        l -> occursin("candidate", get_name(l)),
        collect(get_components(PSY.Line, sys)),
    )

    obj = objective_function(jump_model)
    for line in candidate_lines
        name = get_name(line)
        cost = get(get_ext(line), "project_cost", 0.0)
        JuMP.add_to_expression!(obj, cost, z_var[name])
    end
    set_objective_function(jump_model, obj)
    return
end

"""
    add_candidate_generation_investment_constraints!(decision_model, T) -> variable container

For every component of type `T` whose `ext` dict contains `"is_candidate" => true`:

1. Adds a binary investment variable `x_g ∈ {0,1}` (`GenerationInvestmentVariable`).
2. Adds `ActivePowerVariable[g,t] ≤ p_max_g · x_g` for all time steps
   (`GenerationInvestmentConstraint`).
3. Adds `project_cost · x_g` to the objective (`get_ext(gen)["project_cost"]`,
   defaults to 0 if absent).

Returns the variable container indexed by generator name, or `nothing` if no candidates.
"""
function add_candidate_generation_investment_constraints!(decision_model, T)
    container  = decision_model.internal.container
    jump_model = PSI.get_jump_model(container)
    time_steps = PSI.get_time_steps(container)
    sys        = PSI.get_system(decision_model)

    candidate_gens = filter(
        g -> get(get_ext(g), "is_candidate", false),
        collect(get_components(T, sys)),
    )
    isempty(candidate_gens) && return nothing

    names = get_name.(candidate_gens)

    # Binary investment variable x_g ∈ {0,1}
    x_var = PSI.add_variable_container!(
        container,
        GenerationInvestmentVariable(),
        T,
        names,
    )
    for name in names
        x_var[name] = JuMP.@variable(
            jump_model,
            binary    = true,
            base_name = "GenerationInvestment_$(T)_{$name}",
        )
    end

    # Constraint: p_g[t] ≤ p_max_g * x_g
    p_var = PSI.get_variable(container, ActivePowerVariable(), T)
    con = PSI.add_constraints_container!(
        container,
        GenerationInvestmentConstraint(),
        T,
        names,
        time_steps,
    )
    for gen in candidate_gens
        name  = get_name(gen)
        p_max = PSY.get_max_active_power(gen)
        for t in time_steps
            con[name, t] = JuMP.@constraint(
                jump_model,
                p_var[name, t] <= p_max * x_var[name],
                base_name = "GenerationInvestmentConstraint_$(T)_{$name, $t}",
            )
        end
    end

    # Add project_cost * x_g to the objective
    obj = objective_function(jump_model)
    for gen in candidate_gens
        name = get_name(gen)
        cost = get(get_ext(gen), "project_cost", 0.0)
        JuMP.add_to_expression!(obj, cost, x_var[name])
    end
    set_objective_function(jump_model, obj)

    return x_var
end


"""
    _fc_branch_types(device_models) -> Vector{DataType}

Return the `ACBranch` subtypes present in a device-model dict (the branch types whose
`FlowRateConstraint` bounds receive flow-cancelling shift terms).  DC branches such as
`TwoTerminalGenericHVDCLine` are excluded.
"""
_fc_branch_types(device_models) =
    [T for T in keys(device_models) if T <: PSY.ACBranch]

"""
    build_model_with_flow_canceling_terms(sys; ignore_pf, device_models) -> DecisionModel

Convenience function that builds a complete PTDF model with all flow-cancelling and
investment terms applied in one call.  Internally it:

1. Builds the base PTDF `DecisionModel` via `build!`.
2. Adds binary branch investment variables and flow-cancelling variables for every
   line whose name contains `"candidate"`.
3. Adds big-M linking constraints.
4. Adds PTDF shift terms to existing-line and candidate-line `FlowRateConstraint` bounds,
   for every `ACBranch` type in `device_models` (e.g. `Line`, `TapTransformer`).
5. Adds line investment costs and candidate generation investment constraints to the
   objective.

`device_models` selects the device formulations; pass `RTS_FC_DEVICE_MODELS` for RTS (which
adds `TapTransformer`, renewables, hydro and HVDC). Branch types must use `StaticBranch`.
Existing parallel branches (aggregated `-double_circuit`) are handled automatically.
"""
function build_model_with_flow_canceling_terms(
    sys;
    ignore_pf = true,
    device_models = DEFAULT_FC_DEVICE_MODELS,
    optimizer = DEFAULT_MILP_OPTIMIZER,
)
    ptdf = PTDF(sys)
    name_to_arc_map = _branch_arc_map(ptdf)

    network_model = if ignore_pf
        NetworkModel(PTDFPowerModel; PTDF_matrix = ptdf, use_slacks = true)
    else
        NetworkModel(
            PTDFPowerModel;
            PTDF_matrix = ptdf,
            use_slacks = true,
            power_flow_evaluation = PowerFlows.ACPowerFlow(; calculate_loss_factors = true, calculate_voltage_stability_factors = true),
        )
    end

    template = ProblemTemplate(network_model)
    for (device_type, formulation) in device_models
        set_device_model!(template, device_type, formulation)
    end

    model = DecisionModel(
        template,
        sys;
        optimizer = optimizer,
        name = "UC",
        store_variable_names=true,
    )

    build!(model; output_dir = mktempdir())

    candidate_lines = get_components(x -> contains(x.name, "candidate"), Line, sys)
    z_var = add_branch_investment_variables!(model, candidate_lines, Line)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, Line)
    add_bigM_linking_constraints!(model, candidate_lines, Line, z_var, v_var)

    # Apply shift terms to the existing-branch constraints of every AC branch type.
    for T in _fc_branch_types(device_models)
        add_shift_terms_to_existing_line_constraints!(
            model, T, candidate_lines, v_var, ptdf, sys, name_to_arc_map,
        )
    end

    add_shift_terms_to_candidate_line_constraints!(
        model, sys, candidate_lines, Line, z_var, v_var, ptdf, name_to_arc_map,
    )

    add_candidate_line_investment_costs!(model, z_var)

    add_candidate_generation_investment_constraints!(model, ThermalStandard)

    _fc_add_ptdf_branch_flow_with_fc_expressions!(
        model, sys, ptdf, candidate_lines, v_var, name_to_arc_map,
    )

    return model
end

# ---------------------------------------------------------------------------
# Quadratic loss approximation helpers (self-contained, no voltage scaling)
# ---------------------------------------------------------------------------

"""Variable type for total approximated line losses."""
struct LineLossTotalApproximation <: PSI.VariableType end

"""Constraint type linking bus injections to quadratic losses."""
struct LineLossConstraintApproximation <: PSI.ConstraintType end

"""Expression type for the PTDF branch flow augmented with flow-cancelling correction terms."""
struct PTDFBranchFlowWithFC <: PSI.ExpressionType end

PSI.should_write_resulting_value(::Type{PTDFBranchFlowWithFC}) = true
PSI.convert_result_to_natural_units(::Type{PTDFBranchFlowWithFC}) = true

function _add_existing_branch_fc_expressions!(
    container,
    T::Type,
    time_steps,
    ptdf,
    sys,
    cand_info,
    v_var,
    name_to_arc_map,
)
    key = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlow, T}("")
    if !haskey(container.expressions, key)
        return nothing
    end
    ptdf_branch = container.expressions[key]
    branch_names = axes(ptdf_branch, 1)
    expr = PSI.add_expression_container!(
        container,
        PTDFBranchFlowWithFC(),
        T,
        branch_names,
        time_steps,
    )
    for name in branch_names
        monitored_arc = _resolve_branch_arc(name_to_arc_map, sys, T, name)
        for t in time_steps
            fc_expr = copy(ptdf_branch[name, t])
            for (cand_name, from_num, to_num) in cand_info
                shift = _get_ptdf_shift_term(ptdf, monitored_arc, from_num, to_num)
                JuMP.add_to_expression!(fc_expr, shift, v_var[cand_name, t])
            end
            expr[name, t] = fc_expr
        end
    end
    return expr
end

"""
    _fc_add_ptdf_branch_flow_with_fc_expressions!(model, sys, ptdf, candidate_lines, v_var)

Build PTDFBranchFlowWithFC expression containers for PSY.Line and all existing branch
types (PSY.PhaseShiftingTransformer, PSY.TapTransformer, PSY.Transformer2W) by augmenting
PTDFBranchFlow with FC correction terms. Each entry copies PTDFBranchFlow[name,t] and adds:
  existing arc:   + Σ_cand Δ_{name,cand} · v_cand[t]
  candidate arc:  + (Δ_{name,name}−1) · v_name[t] + Σ_{j≠name} Δ_{name,j} · v_j[t]
Only PSY.Line entries can be candidates; all transformer types always use the existing formula.
A PTDFBranchFlowWithFC container is created only for branch types whose PTDFBranchFlow
expression is present in the container.
"""
function _fc_add_ptdf_branch_flow_with_fc_expressions!(
    model::PSI.DecisionModel,
    sys::PSY.System,
    ptdf,
    candidate_lines,
    v_var,
    name_to_arc_map,
)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)

    candidate_names = Set(get_name.(candidate_lines))
    cand_info = [
        (get_name(c), get_number(get_from(get_arc(c))), get_number(get_to(get_arc(c))))
        for c in candidate_lines
    ]

    # --- PSY.Line expressions ---
    ptdf_line = PSI.get_expression(container, PTDFBranchFlow(), PSY.Line)
    line_names = axes(ptdf_line, 1)

    expr_line = PSI.add_expression_container!(
        container,
        PTDFBranchFlowWithFC(),
        PSY.Line,
        line_names,
        time_steps,
    )

    for name in line_names
        monitored_arc = _resolve_branch_arc(name_to_arc_map, sys, PSY.Line, name)
        is_candidate = name in candidate_names
        for t in time_steps
            expr_line[name, t] = copy(ptdf_line[name, t])

            if is_candidate
                # Candidate arc: self coefficient is (shift - 1), cross coefficient is shift
                for (cand_name, from_num, to_num) in cand_info
                    shift = _get_ptdf_shift_term(ptdf, monitored_arc, from_num, to_num)
                    coeff = cand_name == name ? shift - 1.0 : shift
                    JuMP.add_to_expression!(expr_line[name, t], coeff, v_var[cand_name, t])
                end
            else
                # Existing arc: each candidate contributes shift · v_cand[t]
                for (cand_name, from_num, to_num) in cand_info
                    shift = _get_ptdf_shift_term(ptdf, monitored_arc, from_num, to_num)
                    JuMP.add_to_expression!(expr_line[name, t], shift, v_var[cand_name, t])
                end
            end
        end
    end

    # --- Existing (non-candidate) branch types ---
    for T in (PSY.PhaseShiftingTransformer, PSY.TapTransformer, PSY.Transformer2W, PSY.MonitoredLine)
        _add_existing_branch_fc_expressions!(
            container, T, time_steps, ptdf, sys, cand_info, v_var, name_to_arc_map,
        )
    end

    return expr_line
end

"""
Add continuous loss variables to the model, one per reference bus per time step.
"""
function _fc_add_loss_variables!(model::PSI.DecisionModel)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)
    con_bal = PSI.get_constraint(container, CopperPlateBalanceConstraint(), PSY.System)
    ref_buses = axes(con_bal, 1)

    variable = PSI.add_variable_container!(
        container,
        LineLossTotalApproximation(),
        PSY.System,
        ref_buses,
        time_steps,
    )
    for ref_bus in ref_buses, t in time_steps
        variable[ref_bus, t] = JuMP.@variable(
            PSI.get_jump_model(container),
            base_name = "LineLossTotalApproximation_{$ref_bus}_{$t}"
        )
    end
    return variable
end

"""
Add the loss variable to the copper plate balance (quadratic variant: RHS stays zero).
"""
function _fc_add_loss_to_copperplate_balance!(model::PSI.DecisionModel, loss_variable)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)
    con_bal = PSI.get_constraint(container, CopperPlateBalanceConstraint(), PSY.System)
    ref_buses = axes(con_bal, 1)
    ref_bus = only(ref_buses)
    for t in time_steps
        JuMP.set_normalized_coefficient(con_bal[ref_bus, t], loss_variable[ref_bus, t], 1)
    end
end

"""
    _resistance_by_reduced_name(sys, T, reduced_names, ptdf) -> Dict{String, Float64}

Map each *reduced* branch name (as used on the model axes, possibly `-double_circuit`) of
type `T` to an equivalent resistance.  A direct branch keeps its own `r`; an aggregated
parallel group uses the parallel combination `R_eq = 1 / Σ_i (1/R_i)` over its constituent
PSY components, resolved via `component_to_reduction_name_map`.
"""
function _resistance_by_reduced_name(sys, T, reduced_names, ptdf)
    nrd = ptdf.network_reduction_data
    agg_to_inds = Dict{String, Vector{String}}()
    if haskey(nrd.component_to_reduction_name_map, T)
        for (ind, agg) in nrd.component_to_reduction_name_map[T]
            push!(get!(agg_to_inds, agg, String[]), ind)
        end
    end

    R = Dict{String, Float64}()
    for name in reduced_names
        comp = PSY.get_component(T, sys, name)
        if comp !== nothing
            R[name] = get_r(comp)
        elseif haskey(agg_to_inds, name)
            ginv = sum(1.0 / get_r(PSY.get_component(T, sys, i)) for i in agg_to_inds[name])
            R[name] = 1.0 / ginv
        else
            @warn "No resistance found for reduced branch `$name` ($T); using 0."
            R[name] = 0.0
        end
    end
    return R
end

"""
    _fc_add_quadratic_loss_constraints!(model, sys, ptdf)

Add quadratic loss constraints using the PTDFBranchFlowWithFC expressions:

    loss[ref_bus, t] == -Σ_l R_l · FC_flow[l,t]² - Σ_{non-line branches} R_b · FC_flow[b,t]²

where FC_flow is the PTDFBranchFlowWithFC expression (PTDF flow + FC correction terms).
Covers PSY.Line and all transformer types (PhaseShiftingTransformer, TapTransformer,
Transformer2W) for which a PTDFBranchFlowWithFC expression exists in the container.
Resistances are keyed by reduced axis name (parallel-aggregated where needed).
No voltage scaling is applied — suitable when a flat voltage profile is assumed.
"""
function _fc_add_quadratic_loss_constraints!(model::PSI.DecisionModel, sys::PSY.System, ptdf)
    container = model.internal.container
    time_steps = PSI.get_time_steps(container)

    loss_variable = PSI.get_variable(container, LineLossTotalApproximation(), PSY.System)
    ref_buses = axes(loss_variable, 1)

    constraint = PSI.add_constraints_container!(
        container,
        LineLossConstraintApproximation(),
        PSY.System,
        ref_buses,
        time_steps,
    )

    fc_line = PSI.get_expression(container, PTDFBranchFlowWithFC(), PSY.Line)
    line_names = axes(fc_line, 1)
    R_line = _resistance_by_reduced_name(sys, PSY.Line, line_names, ptdf)

    # Collect (fc_expr, R_dict) for each non-Line branch type with FC expressions
    other_branch_type_data = []
    for T in (PSY.PhaseShiftingTransformer, PSY.TapTransformer, PSY.Transformer2W)
        key = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlowWithFC, T}("")
        if !haskey(container.expressions, key)
            continue
        end
        fc = container.expressions[key]
        R_dict = _resistance_by_reduced_name(sys, T, axes(fc, 1), ptdf)
        push!(other_branch_type_data, (fc, R_dict))
    end

    for ref_bus in ref_buses, t in time_steps
        constraint[ref_bus, t] = JuMP.@constraint(
            PSI.get_jump_model(container),
            loss_variable[ref_bus, t] ==
            -sum(R_line[name] * fc_line[name, t]^2 for name in line_names) -
            sum(
                R_dict[name] * fc[name, t]^2
                for (fc, R_dict) in other_branch_type_data
                for name in axes(fc, 1);
                init = 0.0,
            )
        )
    end
end

"""
    build_model_with_flow_canceling_and_quadratic_losses(sys) -> DecisionModel

Build a PTDF investment model with flow-cancelling constraints **and** a quadratic
PTDF-based transmission loss approximation (no voltage scaling).

This extends `build_model_with_flow_canceling_terms` by adding:
- A loss variable per reference bus per time step.
- The loss variable to the copper plate balance (G - D + loss = 0).
- A quadratic constraint: loss = -Σₖ Rₖ · FC_flow[k,t]²
  where FC_flow[k,t] is the PTDFBranchFlowWithFC expression (PTDF flow + FC correction terms).

Requires an NLP-capable solver (e.g. Ipopt) because of the quadratic constraints.
"""
function build_model_with_flow_canceling_and_quadratic_losses(
    sys;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
    device_models = DEFAULT_FC_DEVICE_MODELS,
)
    ptdf = PTDF(sys)
    name_to_arc_map = _branch_arc_map(ptdf)

    template = ProblemTemplate(NetworkModel(PTDFPowerModel; PTDF_matrix = ptdf, use_slacks = true))
    for (device_type, formulation) in device_models
        set_device_model!(template, device_type, formulation)
    end

    model = DecisionModel(
        template,
        sys;
        optimizer = optimizer,
        name = "UC_QuadLoss",
        store_variable_names = true,
    )

    build!(model; output_dir = mktempdir())

    # --- Flow-cancelling terms (same as build_model_with_flow_canceling_terms) ---
    candidate_lines = get_components(x -> contains(x.name, "candidate"), Line, sys)
    z_var = add_branch_investment_variables!(model, candidate_lines, Line)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, Line)
    add_bigM_linking_constraints!(model, candidate_lines, Line, z_var, v_var)

    for T in _fc_branch_types(device_models)
        add_shift_terms_to_existing_line_constraints!(
            model, T, candidate_lines, v_var, ptdf, sys, name_to_arc_map,
        )
    end
    add_shift_terms_to_candidate_line_constraints!(
        model, sys, candidate_lines, Line, z_var, v_var, ptdf, name_to_arc_map,
    )
    add_candidate_line_investment_costs!(model, z_var)
    add_candidate_generation_investment_constraints!(model, ThermalStandard)

    # --- Quadratic loss approximation ---
    loss_var = _fc_add_loss_variables!(model)
    _fc_add_loss_to_copperplate_balance!(model, loss_var)
    _fc_add_ptdf_branch_flow_with_fc_expressions!(
        model, sys, ptdf, candidate_lines, v_var, name_to_arc_map,
    )
    _fc_add_quadratic_loss_constraints!(model, sys, ptdf)

    return model
end
