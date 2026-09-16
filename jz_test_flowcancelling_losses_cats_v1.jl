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
const PNM = PowerNetworkMatrices

include("SiennaScripts/CircularFlows/mapped_indices.jl")
include("SiennaScripts/CircularFlows/circular_flows.jl")
include("SiennaScripts/add_hvdc.jl")
include("Systems/CATS/build_cats_datacenter.jl")
include("SiennaScripts/build_models.jl")
include("SiennaScripts/build_simulations.jl")
include("SiennaScripts/utils.jl")
include("Systems/5bus/ac_line_expansion_model_example.jl")
include("SiennaScripts/FlowCancelling/build_models.jl")    # flow-cancelling builders
include("SiennaScripts/FlowCancelling/build_models_CATS.jl")    # flow-cancelling builders
include(joinpath(WORKFLOW_DIR, "utils.jl"))
include("Systems/CATS/utils_CATS.jl")
include("Systems/CATS/build_CATS.jl")

const PSY = PowerSystems
const PSI = PowerSimulations

cats_json = joinpath(this_path, "Systems", "CATS", "CATS_saved_reduced_sys.json")

sys = build_cats_system(cats_json)
set_cats_renewable_costs!(sys, 1.0)
scale_cats_loads!(sys, 1.3)
transform_single_time_series!(sys, Hour(2), Hour(2))

output_dir = "./CATS_PTDF"
if !ispath(output_dir)
    mkpath(output_dir)
end

function save_variables_to_csv(results, directory)
    mkpath(directory)
    for (name, values) in PSI.read_variables(results)
        CSV.write(joinpath(directory, "$(name).csv"), values)
    end
end

add_candidate_lines_without_parallel!(sys)
optimizer = optimizer_with_attributes(Gurobi.Optimizer)

system = sys
voltage_limit = 230.0
reduce_radial_branches = true
line_slacks = true

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
set_device_model!(template, DeviceModel(Line, StaticBranch, use_slacks = true, attributes=Dict("filter_function" => x -> get_name(x) in branch_names_voltage_limit)))

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
add_bigM_linking_constraints!(model, candidate_lines, Line, z_var, v_var)

# add flow cancelling to existing line
add_shift_terms_to_existing_line_constraints!(
    model,
    Line,
    candidate_lines,
    v_var,
    ptdf,
    system,
    map_branch_name[Line],
)
# add flow cancelling to candidate lines
add_shift_terms_to_candidate_line_constraints!(
    model,
    system,
    candidate_lines,
    Line,
    z_var,
    v_var,
    ptdf,
    map_branch_name[Line],
)

# add investment costs for candidate lines to objective function
add_candidate_line_investment_costs!(model, z_var)

# remove any slack variables that were added to candidate lines (if line_slacks = true)
if line_slacks
    remove_slack_variables_from_candidate_lines!(model, candidate_lines, Line)
end

loss_var = _fc_add_loss_variables!(model)
_fc_add_loss_to_copperplate_balance!(model, loss_var)
_fc_add_ptdf_branch_flow_with_fc_expressions!(
    model, sys, ptdf, candidate_lines, v_var, map_branch_name[Line],
)
_fc_add_quadratic_loss_constraints!(model, sys, ptdf)

solve!(model)


res_fc_lossless = OptimizationProblemResults(model)
save_variables_to_csv(res_fc_lossless, joinpath(output_dir, "flow_canceling_lossless"))
println("=== Flow-cancelling model solved ===")
println("Objective: ", JuMP.objective_value(model.internal.container.JuMPmodel))

inv_lines_fc = read_variable(res_fc_lossless, PSI.VariableKey{BranchInvestmentVariable, Line}(""))
#inv_gens_fc  = read_variable(res_fc_lossless, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions (with losses):\n", inv_lines_fc)
#println("Generation investment decisions (with losses):\n", inv_gens_fc)

losses_fc = read_variable(res_fc_lossless, PSI.VariableKey{LineLossTotalApproximation, PSY.System}(""))
println("Flow-cancelling losses:\n", losses_fc)