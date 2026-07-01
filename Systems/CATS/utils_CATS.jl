

function remove_slack_variables_from_candidate_lines!(
    decision_model,
    candidate_lines,
    T,
)
    container  = decision_model.internal.container
    time_steps = PSI.get_time_steps(container)
    ub_con = PSI.get_constraint(container, FlowRateConstraint(), T, "ub")
    lb_con = PSI.get_constraint(container, FlowRateConstraint(), T, "lb")

    for line in candidate_lines
        
        k_name   = get_name(line)

        for t in time_steps

            # get constraints
            ub = ub_con[k_name, t]
            lb = lb_con[k_name, t]

            # get slack variables
            slack_ub_var = container.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerSlackUpperBound, Line}("")][k_name,t]
            slack_lb_var = container.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerSlackLowerBound, Line}("")][k_name,t]

            # set coefficients to zero
            set_normalized_coefficient(ub, slack_ub_var, 0.0)
            set_normalized_coefficient(lb, slack_lb_var, 0.0)

            #check
            @assert normalized_coefficient(ub, slack_ub_var) == 0.0
            @assert normalized_coefficient(lb, slack_lb_var) == 0.0

        end
    end
end



function fix_existing_generator_and_source_outputs!(
    model::DecisionModel,
    considered_datetime::DateTime,
    results::OptimizationProblemResults,
    system::System;
    allow_epsilon_adjustment::Bool = false,
    epsilon::Float64 = 1e-6
    )

    printstyled("Fixing existing generator and source outputs to results values for the considered datetime..." * "\n", color=:blue, bold=true)

    # get the base power from results to convert generator outputs from NATURAL_UNITS to SYSTEM_BASE
    results_base_power = results.base_power

    # get the jump model
    container  = model.internal.container

    # get existing generators and sources
    existing_gens = collect(get_components(x -> get(get_ext(x), "candidate", false) == false, Generator, system))
    existing_sources = collect(get_components(x -> get(get_ext(x), "candidate", false) == false, Source, system))

    # cycle through existing generators and sources and fix their outputs to results values for the considered datetime
    for gen_type in unique(typeof.(existing_gens))

        # get subset of generators of type gen_type
        filtered_gens = filter(gen -> typeof(gen) == gen_type, existing_gens)

        # get the results for this gen_type and considered_datetime
        gen_output_results = read_variable(results, PSI.VariableKey{ActivePowerVariable, gen_type}("")) 
        filter!(row -> row[:DateTime] == considered_datetime, gen_output_results)

        # get the variable key and var_array associted with gen_type
        var_key = InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable, gen_type}("")
        var_array = model.internal.container.variables[var_key]

        # cycle through all the generators
        for gen in filtered_gens

            println("Fixing output for generator: ", get_name(gen))

            # get the generator output results for this generator and considered_datetime
            gen_output = filter(row -> row[:name] == get_name(gen), gen_output_results)
            @assert nrow(gen_output) == 1 "Expected exactly one output result for generator $(get_name(gen)) at datetime $(considered_datetime)"

            # convert gen output from NATURAL_UNITS to SYSTEM_BASE
            val_SYSTEM_BASE = gen_output.value[1] / results_base_power
            println("Generator output in SYSTEM_BASE: ", val_SYSTEM_BASE)

            # get the variable for this generator
            var = var_array[get_name(gen), 1]
            
            # fix the variable to the value from results (converted from NATURAL_UNITS to SYSTEM_BASE)
            if allow_epsilon_adjustment
                set_lower_bound(var, val_SYSTEM_BASE - epsilon)
                set_upper_bound(var, val_SYSTEM_BASE + epsilon)
            else
                JuMP.fix(var, val_SYSTEM_BASE; force=true)
                #JuMP.fix(var, 0.0; force=true)
            end

        end
    end

    # Source has one component type but two decision variables: in and out active power.
    if !isempty(existing_sources)

        # get the results for Sources and filter for considered_datetime
        src_in_results = read_variable(results, PSI.VariableKey{ActivePowerInVariable, Source}(""))
        src_out_results = read_variable(results, PSI.VariableKey{ActivePowerOutVariable, Source}(""))
        filter!(row -> row[:DateTime] == considered_datetime, src_in_results)
        filter!(row -> row[:DateTime] == considered_datetime, src_out_results)

        # get the variable keys and var_arrays for Source in and out active power variables
        var_in_key = InfrastructureSystems.Optimization.VariableKey{ActivePowerInVariable, Source}("")
        var_out_key = InfrastructureSystems.Optimization.VariableKey{ActivePowerOutVariable, Source}("")
        var_in_array = model.internal.container.variables[var_in_key]
        var_out_array = model.internal.container.variables[var_out_key]

        # cycle through all existin sources
        for src in existing_sources

            # get Source name
            src_name = get_name(src)

            # find the results associated with src
            src_in_result = filter(row -> row[:name] == src_name, src_in_results)
            src_out_result = filter(row -> row[:name] == src_name, src_out_results)
            @assert nrow(src_in_result) == 1 "Expected exactly one input power result for source $(src_name) at datetime $(considered_datetime)"
            @assert nrow(src_out_result) == 1 "Expected exactly one output power result for source $(src_name) at datetime $(considered_datetime)"

            # convert src in and out power from NATURAL_UNITS to SYSTEM_BASE
            in_val_system_base = src_in_result.value[1] / results_base_power
            out_val_system_base = src_out_result.value[1] / results_base_power

            # fix the variables to the values from results (converted from NATURAL_UNITS to SYSTEM_BASE)
            if allow_epsilon_adjustment
                set_lower_bound(var_in_array[src_name, 1], in_val_system_base - epsilon)
                set_upper_bound(var_in_array[src_name, 1], in_val_system_base + epsilon)
                set_lower_bound(var_out_array[src_name, 1], out_val_system_base - epsilon)
                set_upper_bound(var_out_array[src_name, 1], out_val_system_base + epsilon)
            else
                JuMP.fix(var_in_array[src_name, 1], in_val_system_base; force=true)
                JuMP.fix(var_out_array[src_name, 1], out_val_system_base; force=true)
            end

        end
    end

end


function fix_candidate_generator_outputs!(
    model::DecisionModel,
    considered_datetime::DateTime,
    results::OptimizationProblemResults,
    candidate_gens::Vector{<:Generator};
    allow_epsilon_adjustment::Bool = false,
    epsilon::Float64 = 1e-6,
    constraint::Bool = false
    )

    printstyled("Fixing candidate generator outputs to results values for the considered datetime..." * "\n", color=:blue, bold=true)

    
    # get the base power from results to convert generator outputs from NATURAL_UNITS to SYSTEM_BASE
    results_base_power = results.base_power

    # cycle through existing generators and sources and fix their outputs to results values for the considered datetime
    for gen_type in unique(typeof.(candidate_gens))

        # get subset of generators of type gen_type
        filtered_gens = filter(gen -> typeof(gen) == gen_type, candidate_gens)

        # get the results for this gen_type and considered_datetime
        gen_output_results = read_variable(results, PSI.VariableKey{ActivePowerVariable, gen_type}("")) 
        filter!(row -> row[:DateTime] == considered_datetime, gen_output_results)

        # get the variable key and var_array associted with gen_type
        var_key = InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable, gen_type}("")
        var_array = model.internal.container.variables[var_key]

        # cycle through all the generators
        for gen in filtered_gens

            # get the generator output results for this generator and considered_datetime
            gen_output = filter(row -> row[:name] == get_name(gen), gen_output_results)
            @assert nrow(gen_output) == 1 "Expected exactly one output result for generator $(get_name(gen)) at datetime $(considered_datetime)"

            # convert gen output from NATURAL_UNITS to SYSTEM_BASE
            val_SYSTEM_BASE = gen_output.value[1] / results_base_power

            # get the variable for this generator
            var = var_array[get_name(gen), 1]
            
            if allow_epsilon_adjustment
                set_lower_bound(var, val_SYSTEM_BASE - epsilon)
                set_upper_bound(var, val_SYSTEM_BASE + epsilon)
            elseif constraint
                JuMP.@constraint(jump_model, var == val_SYSTEM_BASE)
            else
                # fix the variable to the value from results (converted from NATURAL_UNITS to SYSTEM_BASE)
                JuMP.fix(var, val_SYSTEM_BASE; force=true)
            end

        end
    end

end


function check_candidate_line_slacks(res)
    
    # get slack variables for the lines
    line_upper_slack = read_variable(res, PSI.VariableKey{FlowActivePowerSlackUpperBound, Line}(""))
    line_lower_slack = read_variable(res, PSI.VariableKey{FlowActivePowerSlackLowerBound, Line}(""))

    # filter to only include lines with nonzero slack values, which indicate a violation of the flow limits
    nonzero_line_upper_slack = filter(:value => f-> abs(f) > 0.0, line_upper_slack)
    nonzero_line_lower_slack = filter(:value => f-> abs(f) > 0.0, line_lower_slack)
    all_line_violations = vcat(nonzero_line_upper_slack, nonzero_line_lower_slack)

    # remove lines that have "intermediate" in the name, which are lines that were added to split existing lines and should not be considered for flow cancelling candidates
    # also filter out any lines with values that are very close to zero
    # filtered_line_violations = filter(row -> abs(row.value) > 1e-6 && !contains(row.name, "intermediate"), all_line_violations)
    filtered_line_violations = filter(row -> abs(row.value) > 1e-6, all_line_violations)
    println("Lines with flow limit violations:")
    println(filtered_line_violations)

    return filtered_line_violations
    
end


function add_new_candidate_line_or_double_existing_line!(
    line_name::String,
    system::System,
    map_branch_name::Dict,
    cost_data::DataFrame,
    num_time_periods::Int
)

    # find the associated arc and existing lines
    arc = find_arc_from_reduced_network(map_branch_name, system, line_name)
    existing_lines_on_arc = collect(get_components(x -> get_arc(x) == arc, Line, system))
    line = first(existing_lines_on_arc)
    @assert line !== nothing "Line $line_name not found in the system."

    if get_ext(line)["split_line"] == true

        println("Line ", get_name(line), " is a split line...")

        # get the associated candidate line
        candidate_line_name = get_ext(line)["associated_candidate_line"]
        candidate_line = get_component(Line, system, candidate_line_name)
        println("Associated candidate line: ", candidate_line_name)

        if get(get_ext(candidate_line), "line_doubled", false) == true
            println("     This line has already been doubled, skipping...")
            return
        end

        # double the cost
        external_field = get_ext(candidate_line)
        println("     Doubling the cost of the associated candidate line, ", candidate_line_name, ": ", external_field["cost"], " -> ", 2*external_field["cost"])
        println("     Doubling the project cost of the associated candidate line, ", candidate_line_name, ": ", external_field["project_cost"], " -> ", 2*external_field["project_cost"])
        external_field["cost"] = 2 * external_field["cost"]  # for AIC code
        external_field["project_cost"] = 2 * external_field["project_cost"]  # for Rodrigo's code, total amoritized cost over the time horizon

        # Double the rating of the line
        println("     Doubling the rating of the line, ", get_name(candidate_line), ": ", get_rating(candidate_line), " -> ", 2*get_rating(candidate_line))
        set_rating!(candidate_line, 2* get_rating(candidate_line))

        # Half the reactance of the line to reflect that it is replacing two parallel lines
        println("     Halving the reactance of the line, ", get_name(candidate_line), ": ", get_x(candidate_line), " -> ", 0.5*get_x(candidate_line))
        set_x!(candidate_line, 1 / (2 ./ get_x(candidate_line)))

        # Mark that this line was doubled
        external_field["line_doubled"] = true 

    else

        # add new candidate
        add_new_candidate_line!(line_name, system, map_branch_name[Line], cost_data, num_time_periods)

    end

end


function get_and_print_candidate_solution_details(res, system)

    ######################################################################
    # Candidate branch details
    ######################################################################

    # read branch investment variable and print number of installed branches
    #inv_l = read_variable(res, InfrastructureSystems.Optimization.VariableKey{BranchInvestmentVariable, Line}(""));
    inv_l = res.variable_values[InfrastructureSystems.Optimization.VariableKey{BranchInvestmentVariable, Line}("")]
    installed_branches = filter(:value => f -> f > 0.5, inv_l)
    installed_branches_names = installed_branches.name

    # print details about the installed branches
    printstyled("\nNumber of branches installed: ", nrow(installed_branches), " out of ", nrow(inv_l), "\n", color=:blue)
    println("Branch investment details:")
    println(installed_branches)

    ######################################################################
    # Candidate generator details
    ######################################################################

    if haskey(res.variable_values, InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable, ThermalStandard}(""))

        # Read generator variables
        thermal_vals = read_variable(res, PSI.VariableKey{ActivePowerVariable, ThermalStandard}(""))

        # Get names of candidate generators
        new_thermal_gen_names = [get_name(gen) for gen in get_components(x -> get(get_ext(x), "candidate", false), ThermalStandard, system)]

        # Filter to only candidate generators
        filtered_thermal_vals = thermal_vals[thermal_vals.name .∈ [new_thermal_gen_names], :]

        # Count thermal generators with zero output for all time periods and with at least one non-zero output
        thermal_zero_count = 0
        thermal_nonzero_count = 0
        built_thermal_gen_names = String[]
        for gen_name in unique(filtered_thermal_vals.name)
            gen_data = filter(:name => n -> n == gen_name, filtered_thermal_vals)
            if all(gen_data.value .== 0.0)
                thermal_zero_count += 1
            else
                thermal_nonzero_count += 1
                push!(built_thermal_gen_names, gen_name)
            end
        end

    else

        new_thermal_gen_names = []
        thermal_zero_count = 0
        thermal_nonzero_count = 0
        built_thermal_gen_names = String[]
        filtered_thermal_vals = DataFrame()

    end


    if haskey(res.variable_values, InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable, RenewableDispatch}("")) 

        # Read generator variables
        renewable_vals = read_variable(res, PSI.VariableKey{ActivePowerVariable, RenewableDispatch}(""))

        # Get names of candidate generators
        new_renewable_gen_names = [get_name(gen) for gen in get_components(x -> get(get_ext(x), "candidate", false), RenewableDispatch, system)]

        # Filter to only candidate generators
        filtered_renewable_vals = renewable_vals[renewable_vals.name .∈ [new_renewable_gen_names], :]

        # Count renewable generators with zero output for all time periods and with at least one non-zero output
        renewable_zero_count = 0
        renewable_nonzero_count = 0
        built_renewable_gen_names = String[]
        for gen_name in unique(filtered_renewable_vals.name)
            gen_data = filter(:name => n -> n == gen_name, filtered_renewable_vals)
            if all(gen_data.value .== 0.0)
                renewable_zero_count += 1
            else
                renewable_nonzero_count += 1
                push!(built_renewable_gen_names, gen_name)
            end
        end

    else

        new_renewable_gen_names = []
        renewable_zero_count = 0
        renewable_nonzero_count = 0
        built_renewable_gen_names = String[]
        filtered_renewable_vals = DataFrame()

    end

    # Build table for candidate generators
    println("\n" * "="^100)
    println("Generator Details:")
    println("="^100)
    
    generator_table = DataFrame(
        Name = String[],
        Type = String[],
        Timestep1 = Float64[],
        MaxActivePower_Timestep1 = Float64[],
        Timestep2 = Float64[],
        MaxActivePower_Timestep2 = Float64[],
        Timestep3 = Float64[],
        MaxActivePower_Timestep3 = Float64[]
    )
    
    # Loop through all candidate generators and populate the table
    for gen_name in vcat(built_thermal_gen_names, built_renewable_gen_names)
        # Get values for first 3 timesteps
        gen_data = filter(:name => n -> n == gen_name, vcat(filtered_thermal_vals, filtered_renewable_vals))
        
        ts1_val = nrow(gen_data) >= 1 ? gen_data[1, :value] : 0.0
        ts2_val = nrow(gen_data) >= 2 ? gen_data[2, :value] : 0.0
        ts3_val = nrow(gen_data) >= 3 ? gen_data[3, :value] : 0.0
        
        # Get max active power from system
        if gen_name in built_thermal_gen_names
            gen = get_component(ThermalStandard, system, gen_name)
            gen_type = "Thermal"
            max_ts1 = get_active_power_limits(gen).max
            max_ts2 = get_active_power_limits(gen).max
            max_ts3 = get_active_power_limits(gen).max
        elseif gen_name in built_renewable_gen_names
            gen = get_component(RenewableDispatch, system, gen_name)
            gen_type = "Renewable"
            max_ts1 = values(get_time_series_array(SingleTimeSeries, gen, "max_active_power"; ignore_scaling_factors=false))[1]
            max_ts2 = values(get_time_series_array(SingleTimeSeries, gen, "max_active_power"; ignore_scaling_factors=false))[2]
            max_ts3 = values(get_time_series_array(SingleTimeSeries, gen, "max_active_power"; ignore_scaling_factors=false))[3]
        else
            error("Generator $gen_name not found in either thermal or renewable candidate lists.")
        end
        
        push!(generator_table, (gen_name, gen_type, ts1_val, max_ts1, ts2_val, max_ts2, ts3_val, max_ts3))
    end
    
    println(generator_table)
    installed_candidate_generator_names = generator_table.Name

    printstyled("Number of candidate generators with zero output for all time periods: ", thermal_nonzero_count + renewable_nonzero_count, " out of ", length(new_thermal_gen_names) + length(new_renewable_gen_names), "\n", color=:blue)
    println("  Thermal: ", thermal_nonzero_count, " out of ", length(new_thermal_gen_names))
    println("  Renewable: ", renewable_nonzero_count, " out of ", length(new_renewable_gen_names))

    return installed_branches_names, installed_candidate_generator_names

end


function remove_lines_not_installed!(
    system::System,
    installed_branches_names::Vector{String}
)

    num_removed_lines = 0

    # Loop through all lines in the system and remove those that are not in the list of installed branche
    for line in get_components(x -> get(get_ext(x), "candidate", false) == true, Line, system)
        if line.name in installed_branches_names
            continue
        else
            remove_component!(system, line)
            num_removed_lines += 1
        end
    end

    printstyled("Removed ", num_removed_lines, " lines that were not installed in the solution.\n", color=:red)

    num_candidate_lines_remaining = length(collect(get_components(x -> get(get_ext(x), "candidate", false) == true, Line, system)))
    printstyled("Number of candidate lines remaining in the system: ", num_candidate_lines_remaining, "\n", color=:blue)

end


function remove_generators_not_installed!(
    system::System,
    installed_candidate_generator_names::Vector{String}
)

    num_generators_removed = 0

    # Loop through all generators in the system and remove those that are not in the list of installed candidate generators
    for gen in get_components(x -> get(get_ext(x), "candidate", false) == true, Generator, system)
        if gen.name in installed_candidate_generator_names
            continue
        else
            remove_component!(system, gen)
            num_generators_removed += 1
        end
    end

    printstyled("Removed ", num_generators_removed, " generators that were not built in the solution.\n", color=:red)

    num_candidate_generators_remaining = length(collect(get_components(x -> get(get_ext(x), "candidate", false) == true, Generator, system)))
    printstyled("Number of candidate generators remaining in the system: ", num_candidate_generators_remaining, "\n", color=:blue)

end


function modify_doubled_lines!(
    system::System,
    orig_system::System,
    installed_branches_names::Vector{String},
)

    for line in get_components(x -> get(get_ext(x), "line_doubled", false) == true, Line, system)

        if get_name(line) in installed_branches_names

            println("Line ", get_name(line), " was doubled to address a flow violation.")

            # get the associated candidate line
            candidate_line = get_component(Line, orig_system, line.name)

            println("Modifying line: ", get_name(candidate_line))
            # modify the line parameters to reflect that this is now a single line instead of two parallel lines

            # double the cost
            external_field = get_ext(candidate_line)
            println("     Doubling the cost of the associated candidate line, ", get_name(candidate_line), ": ", external_field["cost"], " -> ", 2*external_field["cost"])
            println("     Doubling the project cost of the associated candidate line, ", get_name(candidate_line), ": ", external_field["project_cost"], " -> ", 2*external_field["project_cost"])
            external_field["cost"] = 2 * external_field["cost"]  # for AIC code
            external_field["project_cost"] = 2 * external_field["project_cost"]  # for Rodrigo's code, total amoritized cost over the time horizon

            # Double the rating of the line
            println("     Doubling the rating of the line, ", get_name(candidate_line), ": ", get_rating(candidate_line), " -> ", 2*get_rating(candidate_line))
            set_rating!(candidate_line, 2* get_rating(candidate_line))

            # Half the reactance of the line to reflect that it is replacing two parallel lines
            println("     Halving the reactance of the line, ", get_name(candidate_line), ": ", get_x(candidate_line), " -> ", 0.5*get_x(candidate_line))
            set_x!(candidate_line, 1 / (2 ./ get_x(candidate_line)))

            # Mark that this line was doubled
            external_field["line_doubled"] = true 

        end

    end

end


function check_candidate_lines_and_split_lines_not_reduced(
    system::System,
    reduction_data::PNM.NetworkReductionData
)

    for line in get_components(x -> get(get_ext(x), "split_line", false) == true, Branch, system)
        if (line.arc.from.number, line.arc.to.number) in reduction_data.removed_arcs
            error("Line ", get_name(line), " was supposed to be split to avoid reduction, but it was still reduced. Please check the line splitting and reduction logic.")
        end
    end

     for line in get_components(x -> get(get_ext(x), "candidate", false) == true, Line, system)
        if (line.arc.from.number, line.arc.to.number) in reduction_data.removed_arcs
            error("Line ", get_name(line), " was supposed to be split to avoid reduction, but it was still reduced. Please check the line splitting and reduction logic.")
        end
     end

end