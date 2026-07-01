function add_new_line_without_parallel_CATS!(
    sys::System,
    candidate_line::Line
    )

    set_units_base_system!(sys, "SYSTEM_BASE")

    # find the duplicated original line that the candidate line is based on
    # extract the name before the "-i" part
    candidate_line_name = get_name(candidate_line)

    # get all the lines that are parallel to the existing line
    parallel_lines = setdiff(collect(get_components(x -> get_arc(x)==get_arc(candidate_line), Branch, sys)), [candidate_line])
    if length(parallel_lines) == 0
        error("No matching line found for candidate line $(candidate_line_name)")
    end

    # Create intermediate bus: number = from_bus_number + 100000
    intermediate_bus = ACBus(;
        number = length(get_components(ACBus, sys)) + 1, #get_number(bus_from) + 100000,
        name = "intermediate_" * get_arc(candidate_line).from.name * "-" * get_arc(candidate_line).to.name,
        available = true,
        bustype = ACBusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = get_voltage_limits(get_arc(candidate_line).from),
        base_voltage = get_base_voltage(get_arc(candidate_line).from),
        area = get_area(get_arc(candidate_line).from),
        load_zone = get_load_zone(get_arc(candidate_line).from),
    )
    add_component!(sys, intermediate_bus)

    # Check that the base voltage of the intermediate bus matches the from and to buses of the candidate line
    @assert get_base_voltage(intermediate_bus) == get_base_voltage(get_arc(candidate_line).from)
    @assert get_base_voltage(intermediate_bus) == get_base_voltage(get_arc(candidate_line).to)

    # add ext to mark this bus as an intermediate bus for line splitting
    get_ext(intermediate_bus)["added_bus_for_split_arc"] = true

    count = 0
    for existing_line in parallel_lines

        count += 1

        line_name = get_name(existing_line)
        bus_from = get_from(existing_line.arc)
        bus_to   = get_to(existing_line.arc)

        # Build two-segment replacement for the existing line
        seg1 = Line(
            #name = line_name * "_segment_1",
            name = bus_from.name * "-(" * intermediate_bus.name * ")-" * string(count),
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
            name = "(" * intermediate_bus.name * ")-" * bus_to.name * "-" * string(count),
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
        existing_line_rating = get_rating(existing_line)
        remove_component!(sys, existing_line)

        # Add the two-segment replacements
        add_component!(sys, seg1)
        add_component!(sys, seg2)

        # Check that the base voltage of the intermediate bus matches the from and to buses of the candidate line
        @assert get_rating(seg2) == existing_line_rating
        @assert get_rating(seg1) == existing_line_rating

        # add ext to mark these lines as split lines
        get_ext(seg1)["split_line"] = true
        get_ext(seg2)["split_line"] = true

        # associated candidate line
        get_ext(seg1)["associated_candidate_line"] = get_name(candidate_line)
        get_ext(seg2)["associated_candidate_line"] = get_name(candidate_line)

    end

end


function find_arc_from_reduced_network(map_branch_name::Dict, system::System, line_name::String)

    keys_with_line_name = [k for (k, v) in map_branch_name if v == line_name]

    line = get_component(PSY.ACBranch, system, keys_with_line_name[1])

    return get_arc(line)

end


function add_new_candidate_line!(line_name::String, system::System, map_branch_name::Dict, cost_data::DataFrame, num_time_periods::Int)

    # find the associated arc and existing lines
    arc = find_arc_from_reduced_network(map_branch_name, system, line_name)
    existing_lines_on_arc = collect(get_components(x -> get_arc(x) == arc, Line, system))
    reference_line = first(existing_lines_on_arc)

    # set units to be NATURAL_UNUTS
    set_units_base_system!(system, "NATURAL_UNITS")
    
    # name new line
    name = get_name(get_from(arc)) * "-" * get_name(get_to(arc)) * "_candidate"

    # create new line
    new_line = Line(;
        name=name,
        available=true,
        active_power_flow=0.0,
        reactive_power_flow=0.0,
        arc=arc,
        r=0.0,
        x=0.0,
        b=(from=0.0, to=0.0),
        rating=0.0,
        angle_limits=(min=-1.571, max=1.571),
        g=(from=0.0, to=0.0),
        services=Device[],
        ext=Dict{String, Any}(),
    )

    # add new line to system
    add_component!(system, new_line)

    # Set parameters
    set_x!(new_line, 1 / sum((1 ./ get_x(line)) for line in existing_lines_on_arc))
    set_rating!(new_line, sum(get_rating(line) for line in existing_lines_on_arc))

    println("x: ", get_x(new_line))
    println("rating: ", get_rating(new_line))

    # Add external fields
    external_field = get_ext(new_line)
    external_field["candidate"] = true # for AIC code
    external_field["cost"] = cost_data[findfirst(cost_data.line_name .== get_name(reference_line)), :line_cost_millions]*10^6  # for AIC code
    external_field["is_candidate"] = true  # for Rodrigo's code
    external_field["project_cost"] = cost_data[findfirst(cost_data.line_name .== get_name(reference_line)), :amoritized_cost_per_hour] * num_time_periods # for Rodrigo's code, total amoritized cost over the time horizon

    # Split the parallel line by adding a new bus in between
    add_new_line_without_parallel_CATS!(system, new_line)

end


"""
Specialized version of `build_model_with_flow_canceling_terms` for the CATS system
"""
function build_model_with_flow_canceling_terms_CATS(
    system::PSY.System,
    optimizer::MOI.OptimizerWithAttributes,
    num_time_periods::Int;
    voltage_limit::Float64 = 0.0,
    reduce_radial_branches::Bool = false,
    line_slacks::Bool = false,
    candidate_slacks::Bool = false,
    transformer_slacks::Bool = false,
    global_slacks::Bool = false,
    Big_M::Float64 = M_max
    )

    # calculate PTDF matrix
    if reduce_radial_branches
        ybus = PNM.Ybus(system; network_reductions = PNM.NetworkReduction[PNM.RadialReduction()])
    else
        ybus = PNM.Ybus(system)
    end
    ptdf = PTDF(ybus)

    # Create template
    template = ProblemTemplate(
        PSI.NetworkModel(
            PSI.PTDFPowerModel;
            use_slacks = global_slacks,
            reduce_radial_branches = reduce_radial_branches,
            #PTDF_matrix = ptdf
        ),
    )

    # print the number of buses and lines in the system based on the PTDF axes
    num_buses = length(ptdf.axes[1])
    printstyled("Number of buses in PTDF: $num_buses\n", color = :green, bold = true)
    num_lines = length(ptdf.axes[2])
    printstyled("Number of lines in PTDF: $num_lines\n", color = :green, bold = true)

    # Limit the lines included in the model based on voltage level
    # include only brnaches that have at least one bus with base voltage above the limit
    set_units_base_system!(system, "NATURAL_UNITS")
    arcs_voltage_limit = get_components(x -> get_base_voltage(get_to(x))>= voltage_limit || get_base_voltage(get_from(x))>= voltage_limit, Arc, system)
    branch_names_voltage_limit = get_name.(get_components(x -> get_arc(x) in arcs_voltage_limit, Branch, system))
    #branch_names_voltage_limit = get_name.(get_components(x -> get_arc(x) in arcs_voltage_limit && !get(get_ext(x), "is_candidate", false), Branch, system))
    #candidate_branch_names_voltage_limit = get_name.(get_components(x -> get_arc(x) in arcs_voltage_limit && get(get_ext(x), "is_candidate", false), Branch, system))

    # Set device models for ED
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicDispatch))
    #set_device_model!(template, DeviceModel(ThermalStandard, ThermalNoMinDispatch))
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, HydroReservoir, HydroEnergyModelReservoir)
    set_device_model!(template, DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless))
    set_device_model!(template, RenewableNonDispatch, FixedOutput)
    #set_device_model!(template, DeviceModel(Transformer2W, StaticBranch, use_slacks = transformer_slacks, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    set_device_model!(template, DeviceModel(Source, ImportExportSourceModel, use_slacks = false, attributes=Dict("reservation" => false)))
    set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = line_slacks, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    #set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = line_slacks, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    #set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = candidate_slacks, attributes=Dict("filter_function" => x -> get_name(x) in candidate_branch_names_voltage_limit)))


    # Transform time series
    transform_single_time_series!(
           system,
           Dates.Hour(num_time_periods), # horizon
           Dates.Hour(num_time_periods), # interval
       );

    # Create and build decision model
    model = DecisionModel(
        template,
        system;
        name = "SCED",
        optimizer = optimizer,
        horizon = Dates.Hour(num_time_periods),
        #system_to_file = false,
        #initialize_model = true,
        check_numerical_bounds = false,
        optimizer_solve_log_print = true,
        #direct_mode_optimizer = false,
        #rebuild_model = false,
        store_variable_names = true,
        calculate_conflict = false,
        initial_time = timestamp(get_time_series_array(SingleTimeSeries, collect(get_components(PowerLoad, system))[1], "max_active_power"))[1]
    )
    build!(model; output_dir = mktempdir(; cleanup = true))

    # get mapping from branch name to reduced branch name
    # accounts for parallel, reduced branches 
    # e.g., "bus-8727-bus-6811-i_8502" => "bus-8727-bus-6811-i_8double_circuit"
    reduction_data = PNM.get_network_reduction_data(ybus)
    PNM.populate_branch_maps_by_type!(reduction_data)
    map_branch_name = reduction_data.component_to_reduction_name_map

    # list the arcs that are removed by reduction
    removed_arcs = PNM.get_removed_arcs(reduction_data)

    # get a list of all existing lines and candidate lines
    # filter lines based on voltage limit (e.g., only include those above 100kV)
    # filter lines based on whether they are radial (e.g., only include those that are not radial)
    if reduce_radial_branches
        all_lines = collect(get_components(x -> get_available(x)==true && get_arc(x) in arcs_voltage_limit && !((get_number(get_from(get_arc(x))), get_number(get_to(get_arc(x)))) in removed_arcs), PSY.Line, system))
    else
        all_lines = collect(get_components(x -> get_available(x)==true && get_arc(x) in arcs_voltage_limit, PSY.Line, system))
    end
    candidate_lines = filter(x -> get(get_ext(x), "is_candidate", false), all_lines)
    existing_lines  = setdiff(all_lines, candidate_lines)

    # add new variables
    z_var = add_branch_investment_variables!(model, candidate_lines, Line)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, Line)

    # add big-M constraints linking z and v
    add_bigM_linking_constraints!(model, candidate_lines, Line, z_var, v_var, Big_M)

    # add flow cancelling to existing line
    add_shift_terms_to_existing_line_constraints!(
        model,
        existing_lines, # existing_lines are already filtered based on voltage limit
        Line,
        candidate_lines,
        v_var,
        ptdf,
        system;
        map_branch_name = map_branch_name[Line]
    )

    # add flow cancelling to existing transformers
    #existing_xfrm = get_components(get_available, Transformer2W, system)
    # existing_xfrm = collect(get_components(x -> get_available(x)==true && get_arc(x) in arcs_voltage_limit && !((get_number(get_from(get_arc(x))), get_number(get_to(get_arc(x)))) in removed_arcs), PSY.Transformer2W, system))
    # add_shift_terms_to_existing_line_constraints!(
    #     model,
    #     existing_xfrm,
    #     Transformer2W,
    #     candidate_lines,
    #     v_var,
    #     ptdf,
    #     system;
    #     map_branch_name = map_branch_name[Transformer2W]
    # )

    # add flow cancelling to candidate lines
    add_shift_terms_to_candidate_line_constraints!(
        model,
        system,
        candidate_lines,
        Line,
        z_var,
        v_var,
        ptdf,
    )

    # add investment costs for candidate lines to objective function
    add_candidate_line_investment_costs!(model, z_var)

    # add ivestment costs for candidate generators to objective function
    #add_candidate_generation_investment_constraints!(model, ThermalStandard)

    # remove any slack variables that were added to candidate lines (if line_slacks = true)
    if line_slacks
        remove_slack_variables_from_candidate_lines!(model, candidate_lines, Line)
    end

    return model
end


function build_CATS_power_flow(
    system::PSY.System,
    optimizer::MOI.OptimizerWithAttributes,
    num_time_periods::Int;
    voltage_limit::Float64 = 0.0,
    reduce_radial_branches::Bool = false,
    line_slacks::Bool = false,
    global_slacks::Bool = false,
    )

    # calculate PTDF matrix
    if reduce_radial_branches
        ybus = PNM.Ybus(system; network_reductions = PNM.NetworkReduction[PNM.RadialReduction()])
    else
        ybus = PNM.Ybus(system)
    end
    ptdf = PTDF(ybus)

    # print the number of buses and lines in the system based on the PTDF axes
    num_buses = length(ptdf.axes[1])
    printstyled("Number of buses in PTDF: $num_buses\n", color = :green, bold = true)
    num_lines = length(ptdf.axes[2])
    printstyled("Number of lines in PTDF: $num_lines\n", color = :green, bold = true)

    # Create template
    template = ProblemTemplate(
        PSI.NetworkModel(
            PSI.PTDFPowerModel;
            use_slacks = global_slacks,
            reduce_radial_branches = reduce_radial_branches,
            #PTDF_matrix = ptdf
        ),
    )

    # Limit the lines included in the model based on voltage level
    # include only brnaches that have at least one bus with base voltage above the limit
    set_units_base_system!(system, "NATURAL_UNITS")
    arcs_voltage_limit = get_components(x -> get_base_voltage(get_to(x))>= voltage_limit || get_base_voltage(get_from(x))>= voltage_limit, Arc, system)
    branch_names_voltage_limit = get_name.(get_components(x -> get_arc(x) in arcs_voltage_limit, Branch, system))

    # Set device models for ED
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicDispatch))
    #set_device_model!(template, DeviceModel(ThermalStandard, ThermalNoMinDispatch))
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, HydroReservoir, HydroEnergyModelReservoir)
    set_device_model!(template, DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless))
    set_device_model!(template, RenewableNonDispatch, FixedOutput)
    #set_device_model!(template, DeviceModel(Transformer2W, StaticBranch, use_slacks = transformer_slacks, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    set_device_model!(template, DeviceModel(Source, ImportExportSourceModel, use_slacks = false, attributes=Dict("reservation" => false)))
    set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = line_slacks, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    #set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = line_slacks, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    #set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = candidate_slacks, attributes=Dict("filter_function" => x -> get_name(x) in candidate_branch_names_voltage_limit)))

    # Transform time series
    transform_single_time_series!(
           system,
           Dates.Hour(num_time_periods), # horizon
           Dates.Hour(num_time_periods), # interval
       );

    # Create and build decision model
    model = DecisionModel(
        template,
        system;
        name = "SCED",
        optimizer = optimizer,
        horizon = Dates.Hour(num_time_periods),
        #system_to_file = false,
        #initialize_model = true,
        check_numerical_bounds = false,
        optimizer_solve_log_print = true,
        #direct_mode_optimizer = false,
        #rebuild_model = false,
        store_variable_names = true,
        calculate_conflict = false,
        initial_time = timestamp(get_time_series_array(SingleTimeSeries, collect(get_components(PowerLoad, system))[1], "max_active_power"))[1]
    )
    build!(model; output_dir = mktempdir(; cleanup = true))

    return model
end

