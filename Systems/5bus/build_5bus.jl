function update_generation_costs!(sys)
    for g in get_components(ThermalStandard, sys)
        set_operation_cost!(
            g,
            ThermalGenerationCost(
                CostCurve(LinearCurve(get_proportional_term(g.operation_cost.variable.value_curve))),
                0.0,
                0.0,
                0.0,
            )
        )
    end
end

function update_line_ratings!(sys, thermal_multiplier)
    for i in get_components(Line, sys)
        set_rating!(i, get_rating(i)*thermal_multiplier)
        if i.name=="bus-1-bus-4-i_2" set_rating!(i,1.5) end
        if i.name=="bus-4-bus-10-i_7" set_rating!(i,2.0) end
    end
end


function build_matpower_5bus_with_updated_lines()
    PowerSystemCaseBuilder.clear_all_serialized_systems()
    system_name = "matpower_case5_sys"
    sys = build_system(MatpowerTestSystems, system_name)
    sys_ts = build_system(PSITestSystems, "c_sys5_uc"; add_single_time_series=true)
    load_sys = ["bus2", "bus3", "bus4"]
    load_sys_ts = ["Bus2", "Bus3", "Bus4"]
    for (l, lt) in zip(load_sys, load_sys_ts)
        load = get_component(StaticLoad, sys, l)
        load_ts = get_component(StaticLoad, sys_ts, lt)
        ts_array = get_time_series_array(SingleTimeSeries, load_ts, "max_active_power"; ignore_scaling_factors = true)
        values(ts_array)[1] = 1.0 # Set load to 1.0 p.u. for first time step
        new_ts = SingleTimeSeries(; name = "max_active_power", data = ts_array, scaling_factor_multiplier=get_max_active_power)
        add_time_series!(sys, load, new_ts)
    end
    update_generation_costs!(sys)
    update_line_ratings!(sys, 1.0)
    return sys
end

function candidate_projects_data(sys)
    gen1 = get_component(Generator, sys, "gen-1")
    gen4 = get_component(Generator, sys, "gen-4")
    candidate1_willingness_to_pay = LinearCurve(100.0)
    candidate2_willingness_to_pay = LinearCurve(300.0)
    return [
        ThermalStandard( 
            name = "candidate_thermal_1",
            available = gen1.available,
            status = gen1.status,
            bus = gen1.bus,
            active_power = gen1.active_power,
            reactive_power = gen1.reactive_power,
            rating = gen1.rating,
            active_power_limits = (min = 0.0, max = gen1.active_power_limits.max),
            reactive_power_limits = gen1.reactive_power_limits,
            ramp_limits = gen1.ramp_limits,
            operation_cost = gen1.operation_cost,
            base_power = gen1.base_power,
            time_limits = gen1.time_limits,
            must_run = gen1.must_run,
            prime_mover_type = gen1.prime_mover_type,
            fuel = gen1.fuel,
            time_at_status = gen1.time_at_status,
            ext = Dict(
                "willingness_to_pay" => candidate1_willingness_to_pay,
                "is_candidate" => true,
                "project_cost" => 100.0,
            )
        ),
        ThermalStandard( 
            name = "candidate_thermal_2",
            available = gen4.available,
            status = gen4.status,
            bus = gen4.bus,
            active_power = gen4.active_power,
            reactive_power = gen4.reactive_power,
            rating = gen4.rating,
            active_power_limits = (min = 0.0, max = gen4.active_power_limits.max),
            reactive_power_limits = gen4.reactive_power_limits,
            ramp_limits = gen4.ramp_limits,
            operation_cost = gen4.operation_cost,
            base_power = gen4.base_power,
            time_limits = gen4.time_limits,
            must_run = gen4.must_run,
            prime_mover_type = gen4.prime_mover_type,
            fuel = gen4.fuel,
            time_at_status = gen4.time_at_status,
            ext = Dict(
                "willingness_to_pay" => candidate2_willingness_to_pay,
                "is_candidate" => true,
                "project_cost" => 50000.0,
            )
        ),
    ]
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
        ("bus-4-bus-10-i_7", "candidate_line_1", Dict{String, Any}("is_candidate" => true, "project_cost" => 10000.0)),
        ("bus-1-bus-4-i_2",  "candidate_line_2", Dict{String, Any}("is_candidate" => true, "project_cost" => 3000.0)),
    ]
    for (line_name, candidate_name, candidate_ext) in candidate_info
        existing_line = get_component(Line, sys, line_name)
        add_new_line_without_parallel!(sys, existing_line, candidate_name; candidate_ext = candidate_ext)
    end
end

function reconductoring_candidate_data(sys)
    line_to_modify = get_component(Line, sys, "bus-1-bus-4-i_2")
    modify_lines = [
        Line(
            name = "reconductor_candidate_1",
            available = line_to_modify.available,
            active_power_flow = line_to_modify.active_power_flow,
            reactive_power_flow = line_to_modify.reactive_power_flow,
            arc = line_to_modify.arc,
            r = line_to_modify.r,
            x = 0.75 * line_to_modify.x, # new reactance
            b = line_to_modify.b,
            rating = 2 * line_to_modify.rating, # new rating
            angle_limits = line_to_modify.angle_limits,
            ext = Dict(
                "is_candidate" => true,
                "project_cost" => 10.0,
            ),
        )
    ]
    return []
end