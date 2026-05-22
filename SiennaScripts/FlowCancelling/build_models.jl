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

const DEFAULT_MILP_OPTIMIZER = optimizer_with_attributes(Xpress.Optimizer)
#optimizer_with_attributes(Xpress.Optimizer)

const PSI = PowerSimulations
const PSY = PowerSystems
const M_max = 10.0

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
    _get_ptdf_shift_term(ptdf, sys, line_name, from_bus_num, to_bus_num) -> Float64

Return the shift factor for line `line_name` in `sys` due to injecting at `from_bus_num`
minus the factor due to injecting at `to_bus_num` (eq. 3: Δ_{l,k} = H_{from_k,l} − H_{to_k,l}).
Returns 0.0 with a warning when the branch or a bus is not found in `ptdf`.
"""
function _get_ptdf_shift_term(
    ptdf,
    sys,
    line_name::AbstractString,
    from_bus_num::Int,
    to_bus_num::Int,
)
    bus_lookup = ptdf.lookup[1]   # bus number  -> col index
    br_lookup  = ptdf.lookup[2]   # branch arc  -> row index

    line_arc = get_arc(get_component(PSY.Branch, sys, line_name))
    line_tuple = (get_number(get_from(line_arc)), get_number(get_to(line_arc)))
    if !haskey(br_lookup, line_tuple)
        @warn "Branch with arc '$line_arc' and name `$line_name` not found in PTDF; shift term set to 0."
        return 0.0
    end
    if !haskey(bus_lookup, from_bus_num) || !haskey(bus_lookup, to_bus_num)
        @warn "Bus $from_bus_num or $to_bus_num not found in PTDF; shift term set to 0."
        return 0.0
    end

    line_arc_ix    = br_lookup[line_tuple]
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
        #Big_M    = get_rating(line)
        Big_M = M_max
        for t in time_steps
            constraint_ub[name, t] = JuMP.@constraint(jump_model, v_var[name, t] <= Big_M * (1 - z_var[name]))
            constraint_lb[name, t] = JuMP.@constraint(jump_model, v_var[name, t] >= -Big_M * (1 - z_var[name]))
        end
    end
    return
end

"""
    add_shift_terms_to_existing_line_constraints!(
        decision_model, existing_lines, T, candidate_lines, v_var, ptdf, sys)

For each existing (non-candidate) branch in `existing_lines` (type `T`) and each candidate
line `k`, append `Δ_{l,k} · v_k[t]` to the `FlowRateConstraint` upper and lower bounds
(eq. 3).  `Δ_{l,k} = PTDF[from_k, l] − PTDF[to_k, l]` is computed from `ptdf` and `sys`.
"""
function add_shift_terms_to_existing_line_constraints!(
    decision_model,
    existing_lines,
    T,
    candidate_lines,
    v_var,
    ptdf,
    sys,
)
    container  = decision_model.internal.container
    time_steps = PSI.get_time_steps(container)
    ub_con = PSI.get_constraint(container, FlowRateConstraint(), T, "ub")
    lb_con = PSI.get_constraint(container, FlowRateConstraint(), T, "lb")

    for line in existing_lines
        l_name = get_name(line)
        for candidate in candidate_lines
            k_name   = get_name(candidate)
            arc_k    = get_arc(candidate)
            from_num = get_number(get_from(arc_k))
            to_num   = get_number(get_to(arc_k))

            shift = _get_ptdf_shift_term(ptdf, sys, l_name, from_num, to_num)

            for t in time_steps
                # UB: add  shift · v_k[t]
                set_normalized_coefficient(ub_con[l_name, t], v_var[k_name, t], shift)
                # LB: add shift · v_k[t]
                set_normalized_coefficient(lb_con[l_name, t], v_var[k_name, t], shift)
            end
        end
    end
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

        # Self shift: Δ_{k,k} = H_{k, from_k} - H_{k, to_k}
        self_shift = _get_ptdf_shift_term(ptdf, sys, k_name, k_from, k_to)

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
                cross_shift = _get_ptdf_shift_term(ptdf, sys, k_name, kk_from, kk_to)

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
    build_model_with_flow_canceling_terms(sys) -> DecisionModel

Convenience function that builds a complete PTDF model with all flow-cancelling and
investment terms applied in one call.  Internally it:

1. Builds the base PTDF `DecisionModel` via `build!`.
2. Adds binary branch investment variables and flow-cancelling variables for every
   line whose name contains `"candidate"`.
3. Adds big-M linking constraints.
4. Adds PTDF shift terms to existing-line and candidate-line `FlowRateConstraint` bounds.
5. Adds line investment costs and candidate generation investment constraints to the
   objective.
"""
function build_model_with_flow_canceling_terms(sys; ignore_pf = true)
    ptdf = PTDF(sys)
    all_lines       = collect(get_components(PSY.Line, sys))
    candidate_lines = filter(l -> occursin("candidate", get_name(l)), all_lines)
    existing_lines  = filter(l -> !occursin("candidate", get_name(l)), all_lines)

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
    set_device_model!(template, ThermalStandard, ThermalDispatchNoMin)
    set_device_model!(template, Line, StaticBranch)
    set_device_model!(template, PhaseShiftingTransformer, StaticBranch)
    set_device_model!(template, PowerLoad, StaticPowerLoad)

    model = DecisionModel(
        template,
        sys;
        optimizer = DEFAULT_MILP_OPTIMIZER, #Xpress.Optimizer,
        name = "UC",
        store_variable_names=true,
    )

    build!(model; output_dir = mktempdir())

    candidate_lines = get_components(x -> contains(x.name, "candidate"), Line, sys)
    z_var = add_branch_investment_variables!(model, candidate_lines, Line)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, Line)
    add_bigM_linking_constraints!(model, candidate_lines, Line, z_var, v_var)

    existing_lines = get_components(get_available,Line, sys)
    add_shift_terms_to_existing_line_constraints!(
        model,
        existing_lines,
        Line,
        candidate_lines,
        v_var,
        ptdf,
        sys,
    )
    existing_xfrm = get_components(get_available, PhaseShiftingTransformer, sys)
    add_shift_terms_to_existing_line_constraints!(
        model,
        existing_xfrm,
        PhaseShiftingTransformer,
        candidate_lines,
        v_var,
        ptdf,
        sys,
    )

    add_shift_terms_to_candidate_line_constraints!(
        model,
        sys,
        candidate_lines,
        Line,
        z_var,
        v_var,
        ptdf,
    )

    add_candidate_line_investment_costs!(model, z_var)

    add_candidate_generation_investment_constraints!(model, ThermalStandard)

    _fc_add_ptdf_branch_flow_with_fc_expressions!(model, sys, ptdf, candidate_lines, v_var)

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
    for name in branch_names, t in time_steps
        fc_expr = copy(ptdf_branch[name, t])
        for (cand_name, from_num, to_num) in cand_info
            shift = _get_ptdf_shift_term(ptdf, sys, name, from_num, to_num)
            JuMP.add_to_expression!(fc_expr, shift, v_var[cand_name, t])
        end
        expr[name, t] = fc_expr
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

    for name in line_names, t in time_steps
        expr_line[name, t] = copy(ptdf_line[name, t])

        if name in candidate_names
            # Candidate arc: self coefficient is (shift - 1), cross coefficient is shift
            for (cand_name, from_num, to_num) in cand_info
                shift = _get_ptdf_shift_term(ptdf, sys, name, from_num, to_num)
                if cand_name == name
                    coeff = shift - 1.0
                else
                    coeff = shift
                end
                JuMP.add_to_expression!(expr_line[name, t], coeff, v_var[cand_name, t])
            end
        else
            # Existing arc: each candidate contributes shift · v_cand[t]
            for (cand_name, from_num, to_num) in cand_info
                shift = _get_ptdf_shift_term(ptdf, sys, name, from_num, to_num)
                JuMP.add_to_expression!(expr_line[name, t], shift, v_var[cand_name, t])
            end
        end
    end

    # --- Existing (non-candidate) branch types ---
    for T in (PSY.PhaseShiftingTransformer, PSY.TapTransformer, PSY.Transformer2W, PSY.MonitoredLine)
        _add_existing_branch_fc_expressions!(container, T, time_steps, ptdf, sys, cand_info, v_var)
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
    _fc_add_quadratic_loss_constraints!(model, sys)

Add quadratic loss constraints using the PTDFBranchFlowWithFC expressions:

    loss[ref_bus, t] == -Σ_l R_l · FC_flow[l,t]² - Σ_{non-line branches} R_b · FC_flow[b,t]²

where FC_flow is the PTDFBranchFlowWithFC expression (PTDF flow + FC correction terms).
Covers PSY.Line and all transformer types (PhaseShiftingTransformer, TapTransformer,
Transformer2W) for which a PTDFBranchFlowWithFC expression exists in the container.
No voltage scaling is applied — suitable when a flat voltage profile is assumed.
"""
function _fc_add_quadratic_loss_constraints!(model::PSI.DecisionModel, sys::PSY.System)
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
    R_line = Dict(get_name(l) => get_r(l) for l in get_components(get_available, PSY.Line, sys))

    # Collect (fc_expr, R_dict) for each non-Line branch type with FC expressions
    other_branch_type_data = []
    for T in (PSY.PhaseShiftingTransformer, PSY.TapTransformer, PSY.Transformer2W)
        key = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlowWithFC, T}("")
        if !haskey(container.expressions, key)
            continue
        end
        fc = container.expressions[key]
        R_dict = Dict(get_name(b) => get_r(b) for b in get_components(get_available, T, sys))
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
)
    ptdf = PTDF(sys)

    template = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks = true))
    set_device_model!(template, ThermalStandard, ThermalDispatchNoMin)
    set_device_model!(template, Line, StaticBranch)
    set_device_model!(template, PhaseShiftingTransformer, StaticBranch)
    set_device_model!(template, PowerLoad, StaticPowerLoad)

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

    existing_lines = get_components(get_available, Line, sys)
    add_shift_terms_to_existing_line_constraints!(
        model, existing_lines, Line, candidate_lines, v_var, ptdf, sys,
    )
    existing_xfrm = get_components(get_available, PhaseShiftingTransformer, sys)
    add_shift_terms_to_existing_line_constraints!(
        model, existing_xfrm, PhaseShiftingTransformer, candidate_lines, v_var, ptdf, sys,
    )
    add_shift_terms_to_candidate_line_constraints!(
        model, sys, candidate_lines, Line, z_var, v_var, ptdf,
    )
    add_candidate_line_investment_costs!(model, z_var)
    add_candidate_generation_investment_constraints!(model, ThermalStandard)

    # --- Quadratic loss approximation ---
    loss_var = _fc_add_loss_variables!(model)
    _fc_add_loss_to_copperplate_balance!(model, loss_var)
    _fc_add_ptdf_branch_flow_with_fc_expressions!(model, sys, ptdf, candidate_lines, v_var)
    _fc_add_quadratic_loss_constraints!(model, sys)

    return model
end
