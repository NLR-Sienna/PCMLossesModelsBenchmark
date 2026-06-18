function build_rts_system()
    sys = build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization=true)
    transform_single_time_series!(sys, Hour(1), Hour(1))
    hvdc1 = only(get_components(TwoTerminalGenericHVDCLine, sys))
    set_loss!(hvdc1, LinearCurve(0.0))
    set_reactive_power_limits_from!(hvdc1, (min=-100, max=100))
    set_reactive_power_limits_to!(hvdc1, (min=-100, max=100))
    set_active_power_limits_from!(hvdc1, (min=-100, max=100))
    set_active_power_limits_to!(hvdc1, (min=-100, max=100))
    set_available!(get_component(HydroDispatch, sys, "201_HYDRO_4"), false)
    return sys
end

"""
    _rts_replace_line_with_candidate!(sys, existing_line, candidate_name; candidate_ext)

Intermediate-bus trick (mirrors the 5-bus tooling): add a candidate line on the existing
line's arc, then replace the existing line with two series segments through a new
intermediate bus and remove the original.  This keeps the candidate as the *only* branch on
its arc, so candidate lines are never parallel (no `-double_circuit` candidate).

`existing_line` must be a single-circuit line (no parallel sibling), otherwise the candidate
would end up parallel to the sibling.
"""
function _rts_replace_line_with_candidate!(
    sys,
    existing_line,
    candidate_name;
    candidate_ext = Dict{String, Any}(),
)
    line_name = get_name(existing_line)
    bus_from = get_from(existing_line.arc)
    bus_to   = get_to(existing_line.arc)

    add_component!(sys, Line(
        name = candidate_name,
        available = existing_line.available,
        active_power_flow = existing_line.active_power_flow,
        reactive_power_flow = existing_line.reactive_power_flow,
        arc = existing_line.arc,
        r = existing_line.r,
        x = existing_line.x,
        b = existing_line.b,
        rating = existing_line.rating,
        angle_limits = existing_line.angle_limits,
        ext = candidate_ext,
    ))

    intermediate_bus = ACBus(;
        number = get_number(bus_from) + 100000,
        name = "intermediate_bus_" * line_name,
        available = true,
        bustype = ACBusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = get_voltage_limits(bus_from),
        base_voltage = get_base_voltage(bus_from),
        area = get_area(bus_from),
        load_zone = get_load_zone(bus_from),
    )
    add_component!(sys, intermediate_bus)

    seg1 = Line(
        name = line_name * "_segment_1",
        available = existing_line.available,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = bus_from, to = intermediate_bus),
        r = existing_line.r * 0.5,
        x = existing_line.x * 0.5,
        b = (from = existing_line.b.from, to = 0.0),
        rating = existing_line.rating,
        angle_limits = existing_line.angle_limits,
    )
    seg2 = Line(
        name = line_name * "_segment_2",
        available = existing_line.available,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = intermediate_bus, to = bus_to),
        r = existing_line.r * 0.5,
        x = existing_line.x * 0.5,
        b = (from = 0.0, to = existing_line.b.to),
        rating = existing_line.rating,
        angle_limits = existing_line.angle_limits,
    )

    remove_component!(sys, existing_line)
    add_component!(sys, seg1)
    add_component!(sys, seg2)
    return candidate_name
end

"""
    add_candidate_lines_rts!(sys; n_candidates = 2, voltage_threshold = 100.0,
                             project_costs = Float64[]) -> Vector{String}

Add `n_candidates` candidate transmission lines to the RTS system using the intermediate-bus
trick.  Target lines are selected programmatically: single-circuit (no parallel sibling)
`Line`s whose from-bus base voltage exceeds `voltage_threshold` kV, with distinct from-buses
(so intermediate-bus numbers don't collide), taken in name order for determinism.

Returns the candidate line names.  Each candidate's `ext` carries `is_candidate => true`
and a `project_cost`.
"""
function add_candidate_lines_rts!(
    sys;
    n_candidates::Int = 2,
    voltage_threshold::Float64 = 100.0,
    project_costs::Vector{Float64} = Float64[],
)
    # Identify single-circuit arcs (arc appearing exactly once among Lines).
    arc_count = Dict{Tuple{Int,Int}, Int}()
    for l in get_components(Line, sys)
        a = get_arc(l)
        key = (get_number(get_from(a)), get_number(get_to(a)))
        arc_count[key] = get(arc_count, key, 0) + 1
    end

    candidates_sorted = sort!(collect(get_components(Line, sys)); by = get_name)
    chosen = Line[]
    used_from_buses = Set{Int}()
    for l in candidates_sorted
        a = get_arc(l)
        key = (get_number(get_from(a)), get_number(get_to(a)))
        arc_count[key] == 1 || continue                       # single-circuit only
        get_base_voltage(get_from(a)) > voltage_threshold || continue
        from_num = get_number(get_from(a))
        from_num in used_from_buses && continue               # distinct from-buses
        push!(chosen, l)
        push!(used_from_buses, from_num)
        length(chosen) == n_candidates && break
    end

    length(chosen) == n_candidates || error(
        "RTS: found only $(length(chosen)) eligible single-circuit HV lines for " *
        "candidates (requested $n_candidates).",
    )

    default_cost = 1.0e4
    names = String[]
    for (i, line) in enumerate(chosen)
        cost = i <= length(project_costs) ? project_costs[i] : default_cost
        cand_name = "candidate_line_" * get_name(line)
        _rts_replace_line_with_candidate!(
            sys, line, cand_name;
            candidate_ext = Dict{String, Any}("is_candidate" => true, "project_cost" => cost),
        )
        push!(names, cand_name)
    end
    return names
end

function set_renewable_costs!(sys, cost_sign::Float64)
    renewables = collect(get_components(RenewableDispatch, sys))
    for (row, gen) in enumerate(renewables)
        new_cost = RenewableGenerationCost(;
            variable=CostCurve(LinearCurve(cost_sign * 0.01 * row)),
        )
        set_operation_cost!(gen, new_cost)
    end
end

function make_uc_template()
    template = ProblemTemplate(NetworkModel(
        PTDFPowerModel;
        use_slacks=true,
        duals=[CopperPlateBalanceConstraint],
        power_flow_evaluation=PowerFlows.ACPowerFlow(; calculate_loss_factors=true, calculate_voltage_stability_factors = true),
    ))
    set_device_model!(template, Line, StaticBranchBounds)
    set_device_model!(template, TapTransformer, StaticBranchBounds)
    set_device_model!(template, ThermalStandard, ThermalBasicUnitCommitment)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(
        template,
        DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalDispatch),
    )
    return template
end

"""
    build_rts_model_with_quadratic_losses(sys) -> DecisionModel

Build the UC decision model for the RTS system and augment it with the quadratic
P = I²R loss term, without requiring a prior AC power flow solution.

Steps:
1. Build a PTDF UC model using HiGHS (build! only, no solve)
2. Add PSI-managed loss variable + copper-plate balance update
3. Add quadratic constraint: loss = -Σₖ Rₖ (Σⱼ PTDFₖⱼ Pⱼ)²

Solving the resulting nonconvex MIQP requires Gurobi with NonConvex=2.
"""
function build_rts_model_with_quadratic_losses(sys)
    ptdf = PTDF(sys)
    optimizer = optimizer_with_attributes(
        Gurobi.Optimizer,
        "MIPGap" => 0.001,
        "OutputFlag" => 1,
        "DisplayInterval" => 1,   # print every node
        "LogFile" => "gurobi.log",
        "Presolve" => 2,
        "Heuristics" => 0.3,
        "NonConvex" => 2,
        "MIPFocus" => 2,        # bound improvement
        "Cuts" => 2
    )
    model = PSI.DecisionModel(
        make_uc_template(), sys;
        optimizer=optimizer, name="UC", store_variable_names=true,
    )
    PSI.build!(model, output_dir=mktempdir())
    update_copperplate_quadratic_loss_approximation_no_voltage!(model, sys, ptdf)
    return model
end

"""
    build_rts_uc_models_hv(; voltage_threshold = 100.0) -> Dict

Return UC device models for RTS where `Line` and `TapTransformer` carry a
`filter_function` that excludes branches whose **from-bus** base voltage is at
or below `voltage_threshold` kV from the `NetworkFlowConstraint`.

RTS uses `TapTransformer` (not `Transformer2W` as in CATS).
The full system PTDF must still be passed for correct loss accounting.
"""
function build_rts_uc_models_hv(; voltage_threshold::Float64 = 100.0)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    return Dict(
        Line              => DeviceModel(Line, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        TapTransformer    => DeviceModel(TapTransformer, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicUnitCommitment,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_rts_ed_models_hv(; voltage_threshold = 100.0) -> Dict

Return ED device models for RTS with the same LV branch filter as
`build_rts_uc_models_hv`.
"""
function build_rts_ed_models_hv(; voltage_threshold::Float64 = 100.0)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    return Dict(
        Line              => DeviceModel(Line, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        TapTransformer    => DeviceModel(TapTransformer, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_rts_ed_models_acopf() -> Dict

Return ED device models for the RTS system using `ACPPowerModel`.
All branches are modeled with `StaticBranchUnbounded` — no `filter_function` —
because `ACPPowerModel` must see every branch for AC feasibility.
RTS uses `TapTransformer` (not `Transformer2W` as in CATS).
"""
function build_rts_ed_models_acopf()
    return Dict(
        Line                       => StaticBranchUnbounded,
        TapTransformer             => StaticBranchUnbounded,
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end