using Pkg
Pkg.activate(".")

using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using SimpleWeightedGraphs
using Graphs
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using LinearAlgebra
using CSV
using HiGHS
using Ipopt
using Dates
#using Xpress
using Logging
using Gurobi
import PowerSystems as PSY
import PowerSimulations as PSI
import PowerSystemCaseBuilder as PSB

this_path = @__DIR__
kestrel_path = joinpath(this_path, "SiennaScripts", "CircularFlows")
WORKFLOW_DIR = "/projects/wetoowiroc/jzhang2/CATS-CaliforniaTestSystem/Sienna/ak_create_CATS_network_with_candidates"
const PNM = PowerNetworkMatrices

include("SiennaScripts/CircularFlows/mapped_indices.jl")
include("SiennaScripts/CircularFlows/circular_flows.jl")
include("SiennaScripts/add_hvdc.jl")
include("Systems/CATS/build_cats.jl")
include("SiennaScripts/build_models.jl")
include("SiennaScripts/build_simulations.jl")
include("SiennaScripts/utils.jl")
include("Systems/5bus/ac_line_expansion_model_example.jl")
#include("SiennaScripts/FlowCancelling/build_models.jl")    # flow-cancelling builders
include("SiennaScripts/FlowCancelling/build_models_CATS.jl")    # flow-cancelling builders
include(joinpath(WORKFLOW_DIR, "utils.jl"))
include("Systems/CATS/utils_CATS.jl")
include("Systems/CATS/build_CATS.jl")

const PSY = PowerSystems
const PSI = PowerSimulations

optimizer = optimizer_with_attributes(Gurobi.Optimizer)

cats_json = joinpath(this_path, "Systems", "CATS", "CATS_saved_reduced_sys.json")

sys = build_cats_system(cats_json)
set_cats_renewable_costs!(sys, 1.0)
scale_cats_loads!(sys, 2.5)
transform_single_time_series!(sys, Hour(2), Hour(2))

cost_data = CSV.read(joinpath(WORKFLOW_DIR, "CATS_data", "CATS_line_costs_and_lengths.csv"), DataFrame);

system = sys
voltage_limit = 500.0
reduce_radial_branches = true
line_slacks = true
quadratic_losses_voltage_limit = 500.0
max_b = 1 / minimum(get_x.(get_components(Line, system)))
M = 2*pi*max_b

upgraded_arcs, upgraded_arcs_buses, gens_to_duplicate = find_candidate_lines_and_gens(system, 10, 10, 500.0, reduce_radial_branches)

add_candidate_lines_without_parallel!(system, upgraded_arcs,cost_data,1)
add_candidate_generators_to_network!(
           system,
           gens_to_duplicate
       )


model = build_model_with_flow_canceling_and_quadratic_losses_CATS(system;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
    device_models = CAT_FC_MODELS,
    voltage_limit,
    reduce_radial_branches,
    quadratic_losses_voltage_limit,
)

solve!(model)

res_fc_lossless = OptimizationProblemResults(model)
#save_variables_to_csv(res_fc_lossless, joinpath(output_dir, "flow_canceling_lossless"))
println("=== Flow-cancelling model solved ===")
println("Objective: ", JuMP.objective_value(model.internal.container.JuMPmodel))

inv_lines_fc = read_variable(res_fc_lossless, PSI.VariableKey{BranchInvestmentVariable, MonitoredLine}(""))
inv_gens_fc  = read_variable(res_fc_lossless, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions (with losses):\n", inv_lines_fc)
#println("Generation investment decisions (with losses):\n", inv_gens_fc)

losses_fc = read_variable(res_fc_lossless, PSI.VariableKey{LineLossTotalApproximation, PSY.System}(""))
println("Flow-cancelling losses:\n", losses_fc)