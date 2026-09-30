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


function check_if_radial_arc(arc::PSY.Arc, system::PSY.System)

    # get from and to buses of the arc
    from_bus = get_from(arc)
    to_bus = get_to(arc)

    # check if there are any other lines connected to from_bus or to_bus
    from_bus_lines = [x for x in get_components(Arc, system) if get_from(x) == from_bus || get_to(x) == from_bus]
    to_bus_lines = [x for x in get_components(Arc, system) if get_from(x) == to_bus || get_to(x) == to_bus]

    if length(from_bus_lines) == 1 || length(to_bus_lines) == 1
        printstyled("     Arc ", get_name(arc), " is a radial arc.\n", color =:red)
        radial = true
    else
        printstyled("     Arc ", get_name(arc), " is not a radial arc.\n", color=:green)
        radial = false
    end

    return radial

end

function find_line_from_reduced_network(map_branch_name::Dict, system::System, line_name::String)


    keys_with_line_name = [k for (k, v) in map_branch_name if v == line_name]

    if length(keys_with_line_name) == 0
        error("No branch found in the original network corresponding to reduced network line name: ", line_name)
    elseif length(keys_with_line_name) == 1
        line = get_component(PSY.ACBranch, system, keys_with_line_name[1])
    elseif length(keys_with_line_name) > 1
        printstyled("Warning: Multiple branches found in the original network corresponding to reduced network line name ", line_name, "\n", color=:yellow)
        rating = 0.0
        voltage = 0.0
        line = get_component(PSY.ACBranch, system, keys_with_line_name[1])
        for k in keys_with_line_name
            k_line = get_component(PSY.ACBranch, system, k)
            k_rating = get_rating(k_line)
            k_voltage = get_base_voltage(get_from(get_arc(k_line)))
            if k_rating >= rating && k_voltage >= voltage
                rating = k_rating
                voltage = k_voltage
                line = get_component(PSY.ACBranch, system, k)
            end
        end
        printstyled("     Selecting line with highest rating: ", get_name(line), "\n", color=:yellow)
    end

    return line

end

function find_candidate_lines_for_upgrade(
    system::PSY.System,
    results,
    num_lines_to_extract::Int;
    reduce_radial_branches::Bool = true
)

    # get mapping from branch name to reduced network branch name (if network reduction was applied)
    if reduce_radial_branches
        ybus = PNM.Ybus(system; network_reductions = PNM.NetworkReduction[PNM.RadialReduction()])
    else
        ybus = PNM.Ybus(system)
    end
    reduction_data = PNM.get_network_reduction_data(ybus)
    PNM.populate_branch_maps_by_type!(reduction_data)
    map_branch_name = reduction_data.component_to_reduction_name_map

    # # list the arcs that are removed by reduction
    # removed_arcs = PNM.get_removed_arcs(reduction_data)

    # Find the branches with slack variables
    variables = read_variables(results);
    variable_keys = keys(variables);

    # Find non-zero values in upper bound slack variables for lines
    line_upper_slack = variables["FlowActivePowerSlackUpperBound__Line"];
    line_lower_slack = variables["FlowActivePowerSlackLowerBound__Line"];

    # extract only non-zero slack variables for lines
    nonzero_line_upper_slack = filter(:value => f-> abs(f) > 0.0, line_upper_slack)
    nonzero_line_lower_slack = filter(:value => f-> abs(f) > 0.0, line_lower_slack)

    # combine lower and upper slack variables, sort by absolute value
    all_line_violations = vcat(nonzero_line_upper_slack, nonzero_line_lower_slack)
    all_line_violations_sorted = sort(all_line_violations, :value, by=abs, rev=true)

    # Get top 10 unique line names with highest violations
    #upgraded_lines = unique(all_line_violations_sorted.name)[1:num_lines_to_extract]

    # find the arcs corresponding to upgraded lines
    upgraded_lines = []
    upgraded_arcs = []
    upgraded_arcs_buses = []
    i = 0

    while true

        i += 1
        println("i = ", i)

        line_name = all_line_violations_sorted.name[i]

        if ~(line_name in upgraded_lines)

            # find the line with name line_name
            line = find_line_from_reduced_network(map_branch_name[Line], system, line_name)

            # find the arc corresponding to the line
            arc = get_arc(line)

            # check if the arc is radial
            radial = check_if_radial_arc(arc, system)

            println("    Line: ", line_name, " - Arc: ", get_name(arc))
            println("         From bus voltage: ", get_base_voltage(get_from(arc)), " kV - To bus voltage: ", get_base_voltage(get_to(arc)), " kV")
            println("         Is radial arc: ", radial)

            if !radial && arc ∉ upgraded_arcs

                push!(upgraded_lines, line_name)
                push!(upgraded_arcs, arc)
                push!(upgraded_arcs_buses, (get_number(get_from(arc)), get_number(get_to(arc))))
                printstyled("    Number line: ", length(upgraded_lines)+1, "\n", color=:green)

            end

        end

        if length(upgraded_lines) == num_lines_to_extract || i == length(all_line_violations_sorted.name)
            break
        end

    end

    return upgraded_arcs, upgraded_arcs_buses

end

function find_candidate_generators(
    system::System,
    upgraded_arcs::AbstractVector,
    upgraded_arcs_buses::AbstractVector,
    num_buses_to_extract::Int
)

    # Create a vector of length num_lines_to_extract with integer entries
    # The sum should equal num_buses_to_extract
    # Generate random distribution of buses across lines
    buses_per_line = zeros(Int, length(upgraded_arcs))

    # Start by assigning 1 bus to each line (ensuring all lines get at least 0)
    remaining_buses = num_buses_to_extract

    # Randomly distribute the buses
    for i in 1:remaining_buses
        # Choose a random line to assign this bus to
        line_idx = rand(1:length(upgraded_arcs))
        buses_per_line[line_idx] += 1
    end

    # PTDF matrix
    #  PTDF = (A^T × B × A)^(-1) × A^T × B
    PTDF_matrix = PSI.PTDF(system)

    # branch bus numbers
    branch_nums = axes(PTDF_matrix)[2]
    bus_nums = axes(PTDF_matrix)[1]
    branch_nums_indices = [findfirst(x -> x == upgraded_arcs_buses[i], branch_nums) for i in 1:length(upgraded_arcs_buses)]

    # Get PTDF values for the first upgraded line
    bus_nums_top = []

    for i in 1:length(branch_nums_indices)

        buses_per_line[i] > 0 || continue

        println(i)
        println("Branch buses: ", upgraded_arcs_buses[i])

        # find the associated PTDF row for the line
        ptdf_row = PTDF_matrix[branch_nums_indices[i], :]
        sorted_indices = sortperm(abs.(ptdf_row); rev=true)

        # find the first PQ bus with a generator connected to it among the top indices
        count = 0
        added_buses = 0
        for j in sorted_indices[3:end] # skip the top 2 buses since they are usually the from and to buses of the line

            count += 1

            bus = collect(get_components( x -> get_number(x) == bus_nums[j], ACBus, system))[1]

            # keep bus if
            # it is a PQ bus
            # has at least one generator
            # is not already in the list of buses to extract
            if get_bustype(bus) == ACBusTypes.PV && length(get_components(x -> get_bus(x) == bus, Union{PSY.ThermalStandard, PSY.RenewableDispatch}, system)) > 0 && !(bus in bus_nums_top)
                added_buses += 1
                push!(bus_nums_top, bus)
                printstyled("   FOUND, count = ", count, "\n", color=:red)
                printstyled("   j = ", j, "\n", color=:yellow)
                printstyled("   Found bus number: ", bus_nums[j], " - Type: ", get_bustype(bus), " - Connected generators: ", length(get_components(x -> get_bus(x) == bus, Generator, system)), "\n", color=:green)
                if added_buses == buses_per_line[i]
                    break
                end
            end

        end

    end

    # check we extracted the correct number of buses
    length(bus_nums_top) == num_buses_to_extract || error("Number of buses found does not equal number of buses to extract. Found ", length(bus_nums_top), " buses.")

    # inialize generators to duplicate
    gens_to_duplicate = []

    for (i, bus) in enumerate(bus_nums_top)

        all_gens = collect(get_components(x -> get_bus(x) == bus_nums_top[i] && get_available(x) == true, Union{PSY.ThermalStandard, PSY.RenewableDispatch}, system))

        # get a random integer between 1 and length(all_gens)
        # choose random generator to duplicate among the available generators at the bus
        random_index = rand(1:length(all_gens))

        printstyled("Selected generator: ", get_name(all_gens[random_index]), " at bus ", get_number(bus), "\n", color=:green)
        printstyled("     Generator type: ", typeof(all_gens[random_index]), "\n", color=:yellow)

        if all_gens[random_index] in gens_to_duplicate
            printstyled("     Generator already selected for duplication. Skipping.\n", color=:red)

            for j in setdiff(1:length(all_gens), random_index)
                
                if all_gens[j] in gens_to_duplicate
                    printstyled("Selected generator: ", get_name(all_gens[random_index]), " at bus ", get_number(bus), "\n", color=:green)
                    printstyled("     Generator type: ", typeof(all_gens[random_index]), "\n", color=:yellow)
                    printstyled("     Generator already selected for duplication. Skipping.\n", color=:red)
                    continue
                else
                    printstyled("Selected generator: ", get_name(all_gens[random_index]), " at bus ", get_number(bus), "\n", color=:green)
                    printstyled("     Generator type: ", typeof(all_gens[random_index]), "\n", color=:yellow)
                    push!(gens_to_duplicate, all_gens[j])
                    printstyled("     Generator added for duplication.\n", color=:green)
                    break
                end

            end

        else

            push!(gens_to_duplicate, all_gens[random_index])

        end

    end

    # pring the number of generators to duplicate 
    printstyled("Number of generators selected for duplication: ", length(gens_to_duplicate), "\n", color=:blue, bold=true)

    return gens_to_duplicate

end

function find_candidate_lines_and_gens(
        system::PSY.System,
        num_lines_to_extract::Int,
        num_buses_to_extract::Int,
        voltage_limit::Float64 = 0.0,
        reduce_radial_branches::Bool = false
)

    if reduce_radial_branches
        ybus = PNM.Ybus(
            system;
            network_reductions = PNM.NetworkReduction[PNM.RadialReduction()],
        )
    else
        ybus = PNM.Ybus(system)
    end

    ptdf = PTDF(ybus)

    template = ProblemTemplate(
        NetworkModel(
            PTDFPowerModel;
            PTDF_matrix = ptdf,
            use_slacks = true,
            reduce_radial_branches = reduce_radial_branches,
        ),
    )
    device_models = CAT_FC_MODELS
    for (device_type, formulation) in device_models
        set_device_model!(template, device_type, formulation)
    end
    set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = true))

    if voltage_limit != 0
        arcs_voltage_limit = get_components(x -> get_base_voltage(get_to(x))>= voltage_limit || get_base_voltage(get_from(x))>= voltage_limit, Arc, system)
        branch_names_voltage_limit = get_name.(get_components(x -> get_arc(x) in arcs_voltage_limit, Branch, system))
        set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = true, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    end

    optimizer =  optimizer_with_attributes(Gurobi.Optimizer)

    model = DecisionModel(
            template,
            system;
            optimizer = optimizer,
            name = "ED",
            store_variable_names = true,
            calculate_conflict = true,
            optimizer_solve_log_print = true,
        )

    build!(model; output_dir = mktempdir(; cleanup = true))
    solve!(model)
    results = OptimizationProblemResults(model);

    get_optimizer_stats(results)

    upgraded_arcs, upgraded_arcs_buses = find_candidate_lines_for_upgrade(
    system,
    results,
    num_lines_to_extract;
    reduce_radial_branches = reduce_radial_branches,
    );
    gens_to_duplicate =  find_candidate_generators(
        system,
        upgraded_arcs,
        upgraded_arcs_buses,
        num_buses_to_extract
    );
    return upgraded_arcs, upgraded_arcs_buses, gens_to_duplicate
end

function add_candidate_lines_without_parallel!(    
    system::System,
    upgraded_arcs::AbstractVector,
    cost_data::DataFrame,
    num_time_periods::Int
)
    for arc in upgraded_arcs
        name = get_name(get_from(arc)) * "-" * get_name(get_to(arc)) * "_candidate"
        existing_lines_on_arc = collect(get_components(x -> get_name(get_arc(x)) == get_name(arc), Line, system))
        num_lines = length(existing_lines_on_arc)
        if isempty(existing_lines_on_arc)
            printstyled("No existing lines found for arc ", name, ". Skipping candidate line creation.\n", color=:yellow)
            continue
        end
        reference_line = first(existing_lines_on_arc)
        new_line = Line(;
            name=name,
            available=true,
            active_power_flow=0.0,
            reactive_power_flow=0.0,
            arc=reference_line.arc,
            r=reference_line.r/num_lines,
            x=reference_line.x/num_lines,
            b=(from=reference_line.b.from*num_lines, to=reference_line.b.to*num_lines),
            rating=reference_line.rating * num_lines,
            angle_limits=(min=-1.571, max=1.571),
            g=(from=0.0, to=0.0),
            services=Device[],
            ext=Dict{String, Any}(),
        )

        add_component!(system, new_line)
        # Add external fields
        external_field = get_ext(new_line)
        external_field["candidate"] = true # for AIC code
        external_field["cost"] = cost_data[findfirst(cost_data.line_name .== get_name(reference_line)), :line_cost_millions]*10^6 # for AIC code
        external_field["amoritized_cost_per_hour"] = cost_data[findfirst(cost_data.line_name .== get_name(reference_line)), :amoritized_cost_per_hour] # for AIC code
        external_field["is_candidate"] = true  # for Rodrigo's code
        external_field["project_cost"] = cost_data[findfirst(cost_data.line_name .== get_name(reference_line)), :amoritized_cost_per_hour] * num_time_periods # for Rodrigo's code, total amoritized cost over the time horizon

        bus_from = get_from(arc)
        bus_to   = get_to(arc)

        intermediate_bus_area = get_component(Area, system, get_area(bus_from).name)
        intermediate_bus_load_zone = get_component(LoadZone, system, get_load_zone(bus_from).name)
        line_name = bus_from.name * "_" * bus_to.name * "_existing_line"
        existing_bus_numbers = Set(get_number.(get_components(ACBus, system)))
        intermediate_bus_number = get_number(bus_from) + 100000
        while intermediate_bus_number in existing_bus_numbers
            intermediate_bus_number += 100000
        end
        intermediate_bus = ACBus(;
            number = intermediate_bus_number,
            name = "intermediate_bus_" * line_name,
            available = true,
            bustype = ACBusTypes.PQ,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = get_voltage_limits(bus_from),
            base_voltage = get_base_voltage(bus_from),
            area = intermediate_bus_area,
            load_zone = intermediate_bus_load_zone,
        )
        add_component!(system, intermediate_bus)

    # Build two-segment replacement for the existing line
        seg1 = Line(
            name = line_name * "_segment_1",
            available = reference_line.available,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = Arc(; from = bus_from, to = intermediate_bus),
            r = reference_line.r/num_lines * 0.5,
            x = reference_line.x/num_lines * 0.5,
            b = (from = reference_line.b.from*num_lines, to = 0.0),
            rating = reference_line.rating * num_lines,
            angle_limits = reference_line.angle_limits,
        )
        seg2 = Line(
            name = line_name * "_segment_2",
            available = reference_line.available,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = Arc(; from = intermediate_bus, to = bus_to),
            r = reference_line.r/num_lines * 0.5,
            x = reference_line.x/num_lines * 0.5,
            b = (from = 0.0, to = reference_line.b.to*num_lines),
            rating = reference_line.rating * num_lines,
            angle_limits = reference_line.angle_limits,
        )

    # Remove the original line (its arc is still referenced by the candidate line)
        for existing_line in existing_lines_on_arc
            remove_component!(system, existing_line)
        end

    # Add the two-segment replacements
        add_component!(system, seg1)
        add_component!(system, seg2)
    end

end

function add_candidate_generators_to_network!(
    system::System,
    gens_to_duplicate::Vector{Any}
)

    # print the number of gens before adding new gens
    printstyled("Number of gens before adding: ", length(get_components(Generator, system)), "\n", color=:blue, bold=true)   

    #set_units_base_system!(system, "DEVICE_BASE")
    set_units_base_system!(system, "NATURAL_UNITS")

    # Loop through the gens to duplicate and create new gens with the same properties
    for gen in gens_to_duplicate

        name = get_name(gen) * "_candidate"
        println("Adding new generator: ", name, " with max active power: ", get_max_active_power(gen), " MW")
        println("get_rating(gen): ", get_rating(gen))
        gen_bus = get_component(ACBus, system, get_name(gen.bus))

        if typeof(gen) == RenewableDispatch

            new_gen = RenewableDispatch(;
                name=name,
                available=true,
                bus=gen_bus,
                active_power=0.0,
                reactive_power=0.0,
                rating= 0.0,
                prime_mover_type=PrimeMovers.OT,
                reactive_power_limits=nothing,
                power_factor= 0.0,
                operation_cost= RenewableGenerationCost(nothing),
                base_power= get_base_power(gen),
                services=Device[],
                dynamic_injector=nothing,
                ext=Dict{String, Any}(),
            )


        elseif typeof(gen) == ThermalStandard

            new_gen = ThermalStandard(;
                name=name,
                available=true,
                status=get_status(gen),
                bus=gen_bus,
                active_power=0.0,
                reactive_power=0.0,
                rating= 0.0,
                #active_power_limits=(min=get_active_power_limits(gen).min, max=get_active_power_limits(gen).max),
                active_power_limits=(min=0.0, max=0.0),
                reactive_power_limits=nothing,
                ramp_limits=nothing,
                #operation_cost= deepcopy(duplicate_gen_comp.operation_cost), #ThermalGenerationCost(; variable = CostCurve(; value_curve = LinearCurve(linear_operation_cost)), fixed = 0.0, start_up= 0.0, shut_down= 0.0),
                operation_cost= ThermalGenerationCost(nothing),
                base_power= get_base_power(gen),
                ext=Dict{String, Any}(),
            )


        else

            error("Generator type ", typeof(gen), " is not supported in this function. Please update the function to handle this generator type.")

        end
                
        # Add external fields
        external_field = get_ext(new_gen)
        external_field["is_candidate"] = true # for Rodrigo's code, since we are not treating new generators as candidates for this analysis, but this can be updated in the future if needed
        external_field["candidate"] = true # for AIC code
        external_field["project_cost"] = 5000.0

        # Add the new generator to the system
        #println("Number of gens before adding: ", length(get_components(Generator, system)))
        add_component!(system, new_gen)
        #println("Number of gens after adding: ", length(get_components(Generator, system)))

        if typeof(gen) == RenewableDispatch

            println("get_rating(gen): ", get_rating(gen))
            set_rating!(new_gen, get_rating(gen))
            println("get_rating(new_gen): ", get_rating(new_gen))
            set_power_factor!(new_gen, get_power_factor(gen))
            set_active_power!(new_gen, get_max_active_power(gen))
            set_operation_cost!(new_gen, get_operation_cost(gen))
            set_base_power!(new_gen, get_base_power(gen))

            println("gen name: ", get_name(gen))
            println("     rating: ", get_rating(new_gen))
            println("     power factor: ", get_power_factor(new_gen))
            println("     rating * power factor: ", get_rating(new_gen)*get_power_factor(new_gen))

        elseif typeof(gen) == ThermalStandard

            set_active_power!(new_gen, get_active_power(gen))
            set_rating!(new_gen, get_rating(gen))
            set_active_power_limits!(new_gen, (min=0.0, max=get_active_power_limits(gen).max))
            set_operation_cost!(new_gen, get_operation_cost(gen))
            set_base_power!(new_gen, get_base_power(gen))

        else

            error("Generator type ", typeof(gen), " is not supported in this function. Please update the function to handle this generator type.")

        end

        if has_time_series(gen, SingleTimeSeries, "max_active_power")
            # get the time series for the original generator
            gen_time_series = get_time_series_array(SingleTimeSeries, gen, "max_active_power"; ignore_scaling_factors = true)
            
            # Wrap TimeArray in SingleTimeSeries before adding to system
            time_series_data = SingleTimeSeries(
                name = "max_active_power",
                data = gen_time_series;
                scaling_factor_multiplier = get_max_active_power,
            )

            # Add time series to the new generator
            add_time_series!(system, new_gen, time_series_data)

            @assert any(isinf.(values(get_time_series_array(SingleTimeSeries, gen, "max_active_power"; ignore_scaling_factors = false)))) == false
            @assert any(isnan.(values(get_time_series_array(SingleTimeSeries, gen, "max_active_power"; ignore_scaling_factors = true)))) == false
        end

        @assert isinf(get_max_active_power(gen)) == false

    end

    # print the numnber of gens before adding new gens
    printstyled("Number of gens before adding: ", length(get_components(Generator, system)), "\n", color=:green, bold=true)   

end

function build_model_with_flow_canceling_CATS(
    system;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
    device_models = CAT_FC_MODELS,
    voltage_limit::Float64 = 0.0,
    reduce_radial_branches = false,
    quadratic_losses_voltage_limit::Float64 = 0.0,
    )

    candidate_lines_to_convert = collect(get_components(
        x -> occursin("candidate", get_name(x)),
        Line,
        system,
    ))
    for line in candidate_lines_to_convert
        convert_component!(system, line, MonitoredLine)
    end
    # convert_component! replaces each Line in the system, so refresh the collection
    # to hold the newly created MonitoredLine objects.
    candidate_lines = collect(get_components(
        x -> occursin("candidate", get_name(x)),
        MonitoredLine,
        system,
    ))

    # Build the network matrices only after conversion so their reduction maps are
    # keyed by MonitoredLine for the candidate components.
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
            use_slacks = true,
            reduce_radial_branches = reduce_radial_branches,
            PTDF_matrix = ptdf
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

    # Set device models for ED
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicDispatch))
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, HydroReservoir, HydroEnergyModelReservoir)
    set_device_model!(template, DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless))
    set_device_model!(template, RenewableNonDispatch, FixedOutput)
    set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = false, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    set_device_model!(template, DeviceModel(MonitoredLine, StaticBranch, use_slacks = false))

    model = DecisionModel(
        template,
        system;
        name = "SCED",
        optimizer = optimizer,
        check_numerical_bounds = false,
        optimizer_solve_log_print = true,
        store_variable_names = true,
        calculate_conflict = true)
    build!(model; output_dir = mktempdir(; cleanup = true))

    # get mapping from branch name to reduced branch name
    # accounts for parallel, reduced branches 
    # e.g., "bus-8727-bus-6811-i_8502" => "bus-8727-bus-6811-i_8double_circuit"
    reduction_data = PNM.get_network_reduction_data(ybus)
    PNM.populate_branch_maps_by_type!(reduction_data)
    name_to_arc_map = reduction_data.name_to_arc_map

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
    #candidate_lines = filter(x -> get(get_ext(x), "is_candidate", false), all_lines)
    existing_lines  = setdiff(all_lines, candidate_lines)

    # add new variables
    z_var = add_branch_investment_variables!(model, candidate_lines, MonitoredLine)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, MonitoredLine)

    # add big-M constraints linking z and v
    add_bigM_linking_constraints!(model, candidate_lines, MonitoredLine, z_var, v_var)

    # add flow cancelling to existing line
    add_shift_terms_to_existing_line_constraints!(
        model,
        Line,
        candidate_lines,
        v_var,
        ptdf,
        system,
        name_to_arc_map,
    )
    # add flow cancelling to candidate lines
    add_shift_terms_to_candidate_line_constraints!(
        model,
        system,
        candidate_lines,
        MonitoredLine,
        z_var,
        v_var,
        ptdf,
        name_to_arc_map,
    )

    # add investment costs for candidate lines to objective function
    add_candidate_line_investment_costs!(model, z_var)
    add_candidate_generation_investment_constraints!(model, ThermalStandard)
    # # remove any slack variables that were added to candidate lines (if line_slacks = true)
    # if line_slacks
    #     remove_slack_variables_from_candidate_lines!(model, candidate_lines, Line)
    # end

    _fc_add_ptdf_branch_flow_with_fc_expressions!(
        model, system, ptdf, candidate_lines, v_var, name_to_arc_map, quadratic_losses_voltage_limit,
    )
    
    return model
end


function build_model_with_flow_canceling_and_quadratic_losses_CATS(
    system;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
    device_models = CAT_FC_MODELS,
    voltage_limit::Float64 = 0.0,
    reduce_radial_branches = false,
    quadratic_losses_voltage_limit::Float64 = 0.0,
    quadratic_loss_line_names = nothing,
    )

    candidate_lines_to_convert = collect(get_components(
        x -> occursin("candidate", get_name(x)),
        Line,
        system,
    ))
    for line in candidate_lines_to_convert
        convert_component!(system, line, MonitoredLine)
    end
    # convert_component! replaces each Line in the system, so refresh the collection
    # to hold the newly created MonitoredLine objects.
    candidate_lines = collect(get_components(
        x -> occursin("candidate", get_name(x)),
        MonitoredLine,
        system,
    ))

    # Build the network matrices only after conversion so their reduction maps are
    # keyed by MonitoredLine for the candidate components.
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
            use_slacks = true,
            reduce_radial_branches = reduce_radial_branches,
            PTDF_matrix = ptdf
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

    # Set device models for ED
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicDispatch))
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, HydroReservoir, HydroEnergyModelReservoir)
    set_device_model!(template, DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless))
    set_device_model!(template, RenewableNonDispatch, FixedOutput)
    set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = false, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    set_device_model!(template, DeviceModel(MonitoredLine, StaticBranch, use_slacks = false))

    model = DecisionModel(
        template,
        system;
        name = "SCED",
        optimizer = optimizer,
        check_numerical_bounds = false,
        optimizer_solve_log_print = true,
        store_variable_names = true,
        calculate_conflict = true)
    build!(model; output_dir = mktempdir(; cleanup = true))

    # get mapping from branch name to reduced branch name
    # accounts for parallel, reduced branches 
    # e.g., "bus-8727-bus-6811-i_8502" => "bus-8727-bus-6811-i_8double_circuit"
    reduction_data = PNM.get_network_reduction_data(ybus)
    PNM.populate_branch_maps_by_type!(reduction_data)
    name_to_arc_map = reduction_data.name_to_arc_map

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
    #candidate_lines = filter(x -> get(get_ext(x), "is_candidate", false), all_lines)
    existing_lines  = setdiff(all_lines, candidate_lines)

    # add new variables
    z_var = add_branch_investment_variables!(model, candidate_lines, MonitoredLine)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, MonitoredLine)

    # add big-M constraints linking z and v
    add_bigM_linking_constraints!(model, candidate_lines, MonitoredLine, z_var, v_var)

    # add flow cancelling to existing line
    add_shift_terms_to_existing_line_constraints!(
        model,
        Line,
        candidate_lines,
        v_var,
        ptdf,
        system,
        name_to_arc_map,
    )
    # add flow cancelling to candidate lines
    add_shift_terms_to_candidate_line_constraints!(
        model,
        system,
        candidate_lines,
        MonitoredLine,
        z_var,
        v_var,
        ptdf,
        name_to_arc_map,
    )

    # add investment costs for candidate lines to objective function
    add_candidate_line_investment_costs!(model, z_var)

    # # remove any slack variables that were added to candidate lines (if line_slacks = true)
    # if line_slacks
    #     remove_slack_variables_from_candidate_lines!(model, candidate_lines, Line)
    # end

    loss_var = _fc_add_loss_variables!(model)
    _fc_add_loss_to_copperplate_balance!(model, loss_var)
    _fc_add_ptdf_branch_flow_with_fc_expressions!(
        model, system, ptdf, candidate_lines, v_var, name_to_arc_map, quadratic_losses_voltage_limit,
    )
    _fc_add_quadratic_loss_constraints!(
        model,
        system,
        ptdf;
        line_names = quadratic_loss_line_names,
        candidate_line_names = get_name.(candidate_lines),
    )
    
    return model
end

function build_model_with_flow_canceling_and_PWL_quadratic_losses_CATS(
    system;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
    device_models = CAT_FC_MODELS,
    voltage_limit::Float64 = 0.0,
    reduce_radial_branches = false,
    quadratic_losses_voltage_limit::Float64 = 0.0,
    quadratic_loss_line_names = nothing,
    num_pwl_segments::Int = 10,
    )

    candidate_lines_to_convert = collect(get_components(
        x -> occursin("candidate", get_name(x)),
        Line,
        system,
    ))
    for line in candidate_lines_to_convert
        convert_component!(system, line, MonitoredLine)
    end
    # convert_component! replaces each Line in the system, so refresh the collection
    # to hold the newly created MonitoredLine objects.
    candidate_lines = collect(get_components(
        x -> occursin("candidate", get_name(x)),
        MonitoredLine,
        system,
    ))

    # Build the network matrices only after conversion so their reduction maps are
    # keyed by MonitoredLine for the candidate components.
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
            use_slacks = true,
            reduce_radial_branches = reduce_radial_branches,
            PTDF_matrix = ptdf
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

    # Set device models for ED
    set_device_model!(template, DeviceModel(ThermalStandard, ThermalBasicDispatch))
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, DeviceModel(RenewableDispatch, RenewableFullDispatch))
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, HydroReservoir, HydroEnergyModelReservoir)
    set_device_model!(template, DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless))
    set_device_model!(template, RenewableNonDispatch, FixedOutput)
    set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = false, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))
    set_device_model!(template, DeviceModel(MonitoredLine, StaticBranch, use_slacks = false))

    model = DecisionModel(
        template,
        system;
        name = "SCED",
        optimizer = optimizer,
        check_numerical_bounds = false,
        optimizer_solve_log_print = true,
        store_variable_names = true,
        calculate_conflict = true)
    build!(model; output_dir = mktempdir(; cleanup = true))

    # get mapping from branch name to reduced branch name
    # accounts for parallel, reduced branches 
    # e.g., "bus-8727-bus-6811-i_8502" => "bus-8727-bus-6811-i_8double_circuit"
    reduction_data = PNM.get_network_reduction_data(ybus)
    PNM.populate_branch_maps_by_type!(reduction_data)
    name_to_arc_map = reduction_data.name_to_arc_map

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
    #candidate_lines = filter(x -> get(get_ext(x), "is_candidate", false), all_lines)
    existing_lines  = setdiff(all_lines, candidate_lines)

    # add new variables
    z_var = add_branch_investment_variables!(model, candidate_lines, MonitoredLine)
    v_var = add_branch_cancelling_flow_variables!(model, candidate_lines, MonitoredLine)

    # add big-M constraints linking z and v
    add_bigM_linking_constraints!(model, candidate_lines, MonitoredLine, z_var, v_var)

    # add flow cancelling to existing line
    add_shift_terms_to_existing_line_constraints!(
        model,
        Line,
        candidate_lines,
        v_var,
        ptdf,
        system,
        name_to_arc_map,
    )
    # add flow cancelling to candidate lines
    add_shift_terms_to_candidate_line_constraints!(
        model,
        system,
        candidate_lines,
        MonitoredLine,
        z_var,
        v_var,
        ptdf,
        name_to_arc_map,
    )

    # add investment costs for candidate lines to objective function
    add_candidate_line_investment_costs!(model, z_var)

    # # remove any slack variables that were added to candidate lines (if line_slacks = true)
    # if line_slacks
    #     remove_slack_variables_from_candidate_lines!(model, candidate_lines, Line)
    # end

    loss_var = _fc_add_loss_variables!(model)
    _fc_add_loss_to_copperplate_balance!(model, loss_var)
    _fc_add_ptdf_branch_flow_with_fc_expressions!(
        model, system, ptdf, candidate_lines, v_var, name_to_arc_map, quadratic_losses_voltage_limit,
    )
    _fc_add_pwl_loss_constraints!(
        model,
        system,
        ptdf;
        line_names = quadratic_loss_line_names,
        candidate_line_names = get_name.(candidate_lines),
        num_segments = num_pwl_segments,
    )
    
    return model
end
