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

include("Systems/5bus/build_5bus.jl")
include("SiennaScripts/FlowCancelling/build_models.jl")
include("SiennaScripts/utils.jl")
include("Systems/CATS/build_CATS.jl")
include("Systems/CATS/utils_CATS.jl")

# build network
sys = build_matpower_5bus_with_updated_lines()
transform_single_time_series!(sys, Hour(1), Hour(1))
set_available!(get_component(PhaseShiftingTransformer, sys, "bus-3-bus-4-i_5"), false)
candidate_gens = candidate_projects_data(sys)
for gen in candidate_gens
    add_component!(sys, gen)
end
    
# add candidate lines, modidy parallel lines to include an additional bus
add_candidate_line_data_without_parallel!(sys)

# calculate big_M value
set_units_base_system!(sys, "SYSTEM_BASE")
max_b = 1 / minimum(get_x.(get_components(Line, sys)))
Big_M = 2*pi*max_b

# build model with flow cancelling terms for 5-bus system
model = build_model_with_flow_canceling_terms(sys, Big_M; num_time_periods = 1)
jump_model = PSI.get_jump_model(model.internal.container);
write_to_file(jump_model, "model_5-bus" * ".lp");
solve!(model)

# get results
res = OptimizationProblemResults(model)
obj_fun = JuMP.objective_function(model.internal.container.JuMPmodel)

# print solution details
installed_branches_names, installed_candidate_generator_names = get_and_print_candidate_solution_details(res, sys)

# save results and newtork
BASE_DIR = "/Users/akody/Library/CloudStorage/OneDrive-NREL/Projects/2024 LDRD GIQ/PCMLossesModelsBenchmark/Systems/5bus"
# to_json(sys, joinpath(BASE_DIR, "5-bus.json"); force=true);
# serialize_results(res, BASE_DIR)

######################################################################################
# Run power flow for debugging purposes
######################################################################################

# # build the original system
# orig_sys = build_matpower_5bus_with_updated_lines()

# # add generator
# gen1 = get_component(Generator, orig_sys, "gen-1")
# candidate1_willingness_to_pay = LinearCurve(100.0)
# add_component!(orig_sys, ThermalStandard( 
#     name = "candidate_thermal_1",
#     available = gen1.available,
#     status = gen1.status,
#     bus = gen1.bus,
#     active_power = gen1.active_power,
#     reactive_power = gen1.reactive_power,
#     rating = gen1.rating,
#     active_power_limits = (min = 0.0, max = gen1.active_power_limits.max),
#     reactive_power_limits = gen1.reactive_power_limits,
#     ramp_limits = gen1.ramp_limits,
#     operation_cost = gen1.operation_cost,
#     base_power = gen1.base_power,
#     time_limits = gen1.time_limits,
#     must_run = gen1.must_run,
#     prime_mover_type = gen1.prime_mover_type,
#     fuel = gen1.fuel,
#     time_at_status = gen1.time_at_status,
#     ext = Dict(
#         "willingness_to_pay" => candidate1_willingness_to_pay,
#         "is_candidate" => true,
#         "candidate" => true,
#         "project_cost" => 100.0,
#     )
# ))

# # Add candidate line connecting the same buses as the existing line
# existing_line = get_component(Line, orig_sys, "bus-1-bus-4-i_2")
# line_name = get_name(existing_line)
# bus_from = get_from(existing_line.arc)
# bus_to   = get_to(existing_line.arc)
# add_component!(orig_sys, Line(
#     name = "candidate_line_2",
#     available = existing_line.available,
#     active_power_flow = existing_line.active_power_flow,
#     reactive_power_flow = existing_line.reactive_power_flow,
#     arc = existing_line.arc,
#     r = existing_line.r,
#     x = existing_line.x,
#     b = existing_line.b,
#     rating = existing_line.rating,
#     angle_limits = existing_line.angle_limits,
#     ext = Dict(
#         "is_candidate" => true, 
#         "project_cost" => 3000.0
#     )
# ))

# transform_single_time_series!(orig_sys, Hour(1), Hour(1))
# set_available!(get_component(PhaseShiftingTransformer, orig_sys, "bus-3-bus-4-i_5"), false)

# pf_model = build_power_flow_model_5_bus(orig_sys; num_time_periods = 1)

# # fix the dispatches of existing and candidate generators
# # considered_datetime = timestamp(get_time_series_array(SingleTimeSeries, collect(get_components(PowerLoad, orig_sys))[1], "max_active_power"))[1]
# # fix_existing_generator_and_source_outputs!(pf_model, considered_datetime, res, orig_sys; allow_epsilon_adjustment=true, epsilon=1e-1)
# # candidate_gens = collect(get_components(x -> get_name(x) ∈ installed_candidate_generator_names, Generator, orig_sys))
# # fix_candidate_generator_outputs!(pf_model, considered_datetime, res, candidate_gens; allow_epsilon_adjustment=true, epsilon=1e-0)

# # jump_model = PSI.get_jump_model(pf_model.internal.container);
# # write_to_file(jump_model, "pf_model_5-bus" * ".lp");

# solve!(pf_model)