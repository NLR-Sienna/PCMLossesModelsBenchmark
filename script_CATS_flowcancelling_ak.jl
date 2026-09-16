using Pkg
Pkg.activate(".")

using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using Ipopt
using Xpress
using Dates
using JuMP
using CSV
using JLD2
using Logging
using InfrastructureSystems
using TimeSeries
using DataFrames
import PowerNetworkMatrices

const PSI = PowerSimulations
const PSY = PowerSystems
const PNM = PowerNetworkMatrices

include("Systems/CATS/utils_CATS.jl")
include("Systems/CATS/build_CATS.jl")
include("SiennaScripts/FlowCancelling/build_models.jl")
include("SiennaScripts/utils.jl")

# import CATS network with candidate lines, new generators, and new loads
# system_name = "CATS_Sienna_ak_version"
# num_loads = 1000
# select_loads_to_duplicate_randomly = true
# num_lines_to_extract = 22
# gens_to_duplicate = 130
#system = System("./Systems/CATS/CATS_Sienna_ak_version_w_new_lines_and_new_gens.json")
#system = System("./Systems/CATS/CATS_Sienna_ak_" * string(num_loads) * (select_loads_to_duplicate_randomly ? "_rand" : "") * "_loads_" * string(num_lines_to_extract) * "_lines_" * string(gens_to_duplicate) * "_gens.json")
# system = System("./Systems/CATS/summer_peak_only/CATS_Sienna_ak_" * string(num_loads) * (select_loads_to_duplicate_randomly ? "_rand" : "") * "_loads_" * string(num_lines_to_extract) * "_lines_" * string(gens_to_duplicate) * "_gens_SP_only.json")
#system = System("./Systems/CATS/imported/summer_peak_only/CATS_Sienna_ak_SP_only_2026-05-28_10-18-58.json")
date = "2026-06-09_14-07-09"
system = System("./Systems/CATS/imported/summer_peak_only/CATS_Sienna_ak_SP_only_" * date * ".json")

# load line cost data
cost_data = CSV.read("./Systems/CATS/data/CATS_line_costs_and_lengths.csv", DataFrame)

# extract the number of time periods, so we can calculate the total amoritized cost over the time horizon for each line
#num_time_periods = length(timestamp(get_time_series_array(SingleTimeSeries, collect(get_components(PowerLoad, system))[1], "max_active_power")))
num_time_periods = 1

# Each candidate line is parallel to an existing line
# Split existing line into two segments to avoid PTDF reduction issues, where lines would be merged
for candidate_line in get_components(x -> get(get_ext(x), "is_candidate", false), Line, system)
    println("Adding candidate line: ", get_name(candidate_line))
    add_new_line_without_parallel_CATS!(system, candidate_line)
end

# select optimizer with memory/performance settings
optimizer = optimizer_with_attributes(
    Xpress.Optimizer,
    "MIPRELSTOP" => 0.0,
    #"OUTPUTLOG" => 0,  # Reduce output to save memory
)

# reduce the radial branches in the simulation
reduce_radial_branches = true

# set voltage limit to include lines above this limit
voltage_limit = 115.0

# number of time periods
num_time_periods = 1

# initialize
count = 0
res = nothing
model = nothing
installed_branches_names = nothing 
installed_candidate_generator_names = nothing

# while loop to add new line capacity if violations are found
while true

    count += 1
    println("\nIteration ", count, ": Solving model with flow cancelling terms...")

    # set units to be NATURAL_UNUTS
    #set_units_base_system!(system, "NATURAL_UNITS")

    # calculate big_M value
    set_units_base_system!(system, "SYSTEM_BASE")
    max_b = 1 / minimum(get_x.(get_components(Line, system)))
    M = 2*pi*max_b

    # build model with flow cancelling terms for CATS system
    model = build_model_with_flow_canceling_terms_CATS(
        system,
        optimizer,
        num_time_periods;
        voltage_limit = voltage_limit,
        reduce_radial_branches = reduce_radial_branches,
        line_slacks = true,
        candidate_slacks = false,
        transformer_slacks = true,
        global_slacks = false,
        Big_M = M
        )
    # jump_model = PSI.get_jump_model(pf_model.internal.container);
    # write_to_file(jump_model, "pf_model_" * system_name * "_" * string(considered_datetime) * ".lp");
    solve!(model)

    # get results
    res = OptimizationProblemResults(model)

    # determine which lines have flow limit violations based on the slack variable values
    filtered_line_violations = check_candidate_line_slacks(res)

    # calculate the map from orginal line names to reduced line names
    if reduce_radial_branches
        ybus = PNM.Ybus(system; network_reductions = PNM.NetworkReduction[PNM.RadialReduction()])
    else
        ybus = PNM.Ybus(system)
    end
    reduction_data = PNM.get_network_reduction_data(ybus)
    PNM.populate_branch_maps_by_type!(reduction_data)
    map_branch_name = reduction_data.component_to_reduction_name_map

    # Check that candidate lines and split lines were not reduced
    check_candidate_lines_and_split_lines_not_reduced(system, reduction_data)

   # Add new candidate branches
    for line_name in filtered_line_violations.name
        printstyled("Line with violation: ", line_name; color=:red)
        add_new_candidate_line_or_double_existing_line!(line_name, system, map_branch_name[Line], cost_data, num_time_periods)
    end

    # if there are no violations, break from loop and print details about the candidate solution
    if nrow(filtered_line_violations) == 0
        printstyled("No more flow limit violations, stopping iteration."; color=:red)
        installed_branches_names, installed_candidate_generator_names = get_and_print_candidate_solution_details(res, system)
        break
    end

end

# save new system with split lines
# BASE_DIR = "/Users/akody/Library/CloudStorage/OneDrive-NREL/Projects/2024 LDRD GIQ/PCMLossesModelsBenchmark/Systems/CATS/results/summer_peak_only_with_split_lines/" * string(date)
# to_json(system, joinpath(BASE_DIR, "CATS_Sienna_ak_SP_only_split_lines" * "_" * date * ".json"); force=true);

# save results
#PSI.serialize_results(res, "/Users/akody/Library/CloudStorage/OneDrive-NREL/Projects/2024 LDRD GIQ/PCMLossesModelsBenchmark/Systems/CATS/results/summer_peak_only_with_split_lines/" * string(date))

######################################################################################
# Run power flow to test solution
######################################################################################

# get original system before new buses / split lines were added
orig_system = System("./Systems/CATS/imported/summer_peak_only/CATS_Sienna_ak_SP_only_" * date * ".json")

# remove lines that were not selected for installation and candidate generators that were not selected for installation from the system
remove_lines_not_installed!(orig_system, installed_branches_names)
remove_generators_not_installed!(orig_system, installed_candidate_generator_names)

# split the lines that are parallel to the selected candidate lines
# for candidate_line in get_components(x -> get(get_ext(x), "is_candidate", false), Line, orig_system)
#     println("Adding candidate line: ", get_name(candidate_line))
#     add_new_line_without_parallel_CATS!(orig_system, candidate_line)
# end

# double the lines that were doubled due to power flow violations
modify_doubled_lines!(system, orig_system, installed_branches_names)

# generate power flow model
pf_model = build_CATS_power_flow(
    orig_system,
    optimizer,
    1;
    voltage_limit = 115.0,
    reduce_radial_branches = true,
    line_slacks = false,
    global_slacks = false,
    )

# fix the dispatches of existing and candidate generators
considered_datetime = timestamp(get_time_series_array(SingleTimeSeries, collect(get_components(PowerLoad, orig_system))[1], "max_active_power"))[1]
fix_existing_generator_and_source_outputs!(pf_model, considered_datetime, res, orig_system)
candidate_gens = collect(get_components(x -> get_name(x) ∈ installed_candidate_generator_names, Generator, orig_system))
fix_candidate_generator_outputs!(pf_model, considered_datetime, res, candidate_gens)

jump_model = PSI.get_jump_model(pf_model.internal.container);
write_to_file(jump_model, "pf_model_" * "CATS" * "_" * string(considered_datetime) * ".lp");

solve!(pf_model)

# get results
res = OptimizationProblemResults(pf_model)

# save new system with candidate lines added
# BASE_DIR = "/Users/akody/Library/CloudStorage/OneDrive-NREL/Projects/2024 LDRD GIQ/PCMLossesModelsBenchmark/Systems/CATS/results/summer_peak_only_with_upgrades_and_generators/" * string(date)
# to_json(system, joinpath(BASE_DIR, "CATS_Sienna_ak_with_upgrades_and_generators" * "_" * date * ".json"); force=true);

# save results
# PSI.serialize_results(res, "/Users/akody/Library/CloudStorage/OneDrive-NREL/Projects/2024 LDRD GIQ/PCMLossesModelsBenchmark/Systems/CATS/results/summer_peak_only_with_upgrades_and_generators/" * string(date))