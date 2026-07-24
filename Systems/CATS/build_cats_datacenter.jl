# CATS-specific device models.
# CATS uses Transformer2W (not TapTransformer) and has SynchronousCondenser.
const CATS_UC_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicUnitCommitment,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

const CATS_ED_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicDispatch,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    SynchronousCondenser => SynchronousCondenserBasicDispatch,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

"""
    build_cats_uc_models_hv(; voltage_threshold = 100.0) -> Dict

Return UC device models for CATS where `Line` and `Transformer2W` carry a
`filter_function` that excludes branches whose **from-bus** base voltage is at
or below `voltage_threshold` kV from the `NetworkFlowConstraint`.

The full system PTDF must still be passed to the PSI `NetworkModel`; this filter
only removes the per-branch flow-limit constraints for LV branches — it does NOT
remove their contribution to the quadratic loss term.
"""
function build_cats_uc_models_hv(; voltage_threshold::Float64 = 100.0, bounded = false)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    if bounded
        branch_model = StaticBranchBounds
    else
        branch_model = StaticBranchUnbounded
    end
    return Dict(
        Line          => DeviceModel(Line, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        Transformer2W => DeviceModel(Transformer2W, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicUnitCommitment,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_cats_ed_models_hv(; voltage_threshold = 100.0) -> Dict

Return ED device models for CATS with the same LV branch filter as
`build_cats_uc_models_hv`. Includes `SynchronousCondenser` which is
present only in the CATS ED model.
"""
function build_cats_ed_models_hv(; voltage_threshold::Float64 = 100.0, bounded = false)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    if bounded
        branch_model = StaticBranchBounds
    else
        branch_model = StaticBranchUnbounded
    end
    return Dict(
        Line          => DeviceModel(Line, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        Transformer2W => DeviceModel(Transformer2W, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        SynchronousCondenser       => SynchronousCondenserBasicDispatch,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_cats_system(cats_json_path) -> System

Load the CATS system from the saved JSON, add two synthetic HVDC links
(Newark–NRS and Metcalf–SanJoseB at ±1000 MW, zero loss), and transform
the time series to 1-hour horizon / 1-hour interval.

The HVDC lines are required to create inter-area shortcuts that the
optimizer can exploit, making circular flows possible.
"""
function build_cats_system(cats_json_path::String)
    sys = System(cats_json_path; runchecks = false)
    transform_single_time_series!(sys, Hour(1), Hour(1))
    add_internal_hvdc!(sys)
    return sys
end

"""
    set_cats_renewable_costs!(sys, cost_sign)

Scale renewable generator costs.
cost_sign > 0 → small positive cost (HVDC loop still profitable → circular flow scenario).
cost_sign < 0 → small negative cost (over-generation incentivized, no HVDC arbitrage → baseline).
"""
function set_cats_renewable_costs!(sys::PSY.System, cost_sign::Float64)
    renewables = collect(PSY.get_components(RenewableDispatch, sys))
    for (row, gen) in enumerate(renewables)
        new_cost = RenewableGenerationCost(;
            variable = CostCurve(LinearCurve(cost_sign * 0.01 * row)),
        )
        PSY.set_operation_cost!(gen, new_cost)
    end
end

function scale_cats_loads!(sys::PSY.System, scale_factor::Float64)
    loads = collect(PSY.get_components(PowerLoad, sys))
    for load in loads
        max_p = PSY.get_max_active_power(load)
        p = PSY.get_active_power(load)
        PSY.set_active_power!(load, p * scale_factor)
        PSY.set_max_active_power!(load, max_p * scale_factor)
    end
end

"""
    build_cats_ed_models_acopf() -> Dict

Return ED device models for the CATS system using `ACPPowerModel`.
All branches are modeled with `StaticBranchUnbounded` — no `filter_function` —
because `ACPPowerModel` must see every branch for AC feasibility.
Losses are implicit in the nonlinear AC formulation; no separate loss term is needed.
"""
function build_cats_ed_models_acopf(;bounded = false)
    if bounded
        branch_model = StaticBranchBounds
    else
        branch_model = StaticBranchUnbounded
    end
    return Dict(
        Line                       => branch_model,
        Transformer2W              => branch_model,
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        SynchronousCondenser       => SynchronousCondenserBasicDispatch,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

function add_candidate_line_data!(sys)
    line1 = get_component(Line, sys, "bus-4-bus-10-i_7")
    line2 = get_component(Line, sys, "bus-1-bus-4-i_2")
    new_candidates = [
        Line(
            name = "candidate_line_1",
            available = line1.available,
            active_power_flow = line1.active_power_flow,
            reactive_power_flow = line1.reactive_power_flow,
            arc = line1.arc,
            r = line1.r,
            x = line1.x,
            b = line1.b,
            rating = line1.rating,
            angle_limits = line1.angle_limits,
            ext = Dict(
                "is_candidate" => true,
                "project_cost" => 100000.0,
            ),
        ),
        Line(
            name = "candidate_line_2",
            available = line2.available,
            active_power_flow = line2.active_power_flow,
            reactive_power_flow = line2.reactive_power_flow,
            arc = line2.arc,
            r = line2.r,
            x = line2.x,
            b = line2.b,
            rating = line2.rating,
            angle_limits = line2.angle_limits,
            ext = Dict(
                "is_candidate" => true,
                "project_cost" => 30000.0,
            ),
        ),
    ]
    for line in new_candidates
        add_component!(sys, line)
    end

end

function add_new_line_without_parallel!(sys, existing_line, candidate_name; candidate_ext = Dict{String, Any}())
    line_name = get_name(existing_line)
    bus_from = get_from(existing_line.arc)
    bus_to   = get_to(existing_line.arc)

    # Add candidate line connecting the same buses as the existing line
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

    # Create intermediate bus: number = from_bus_number + 100000
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

    # Build two-segment replacement for the existing line
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

    # Remove the original line (its arc is still referenced by the candidate line)
    remove_component!(sys, existing_line)

    # Add the two-segment replacements
    add_component!(sys, seg1)
    add_component!(sys, seg2)
end

function update_line_to_parallel!(sys, existing_line)
    line_name = PowerSystems.get_name(existing_line)
    bus_from = get_from(existing_line.arc)
    bus_to   = get_to(existing_line.arc)

    # Create intermediate bus: number = from_bus_number + 100000
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

    # Build two-segment replacement for the existing line
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

    # Add the two-segment replacements
    add_component!(sys, seg1)
    add_component!(sys, seg2)

    # Remove the original line (its arc is still referenced by the candidate line)
    remove_component!(sys, existing_line)
end


function add_candidate_line_data_without_parallel!(sys)
    candidate_info = [
        ("bus-1951-bus-2835-i_980", "candidate_line_1", Dict{String, Any}("is_candidate" => true, "project_cost" => 10000.0)),
        ("bus-2835-bus-7636-i_1912",  "candidate_line_2", Dict{String, Any}("is_candidate" => true, "project_cost" => 3000.0)),
    ]
    for (line_name, candidate_name, candidate_ext) in candidate_info
        existing_line = get_component(Line, sys, line_name)
        add_new_line_without_parallel!(sys, existing_line, candidate_name; candidate_ext = candidate_ext)
    end
end
