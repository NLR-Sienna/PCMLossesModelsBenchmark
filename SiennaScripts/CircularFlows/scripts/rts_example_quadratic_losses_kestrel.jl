# NOTE: Requires Gurobi with NonConvex=2.
# Deploy at project root on Kestrel, then run:
#   julia rts_example_quadratic_losses_kestrel.jl
using Pkg
this_path = @__DIR__
Pkg.activate(this_path)
Pkg.instantiate()
kestrel_path = joinpath(this_path, "SiennaScripts", "CircularFlows")

using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using HiGHS
using Gurobi
using Graphs
using SimpleWeightedGraphs
using Dates
using DataFrames

include(joinpath(kestrel_path, "mapped_indices.jl"))
include(joinpath(kestrel_path, "circular_flows.jl"))
include(joinpath(kestrel_path, "../../Systems/RTS/build_rts.jl"))
include(joinpath(kestrel_path, "../build_models.jl"))
include(joinpath(kestrel_path, "../utils.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

println("=== Scenario A: with circular flow (positive RE cost_sign) ===")
sys_a = build_rts_system()
set_renewable_costs!(sys_a, 1.0)
model_a = build_rts_model_with_quadratic_losses(sys_a)
PSI.solve!(model_a)
results_a = PSI.OptimizationProblemResults(model_a)
res_vars_a = results_a.variable_values
data_a = PSI.get_power_flow_data(
    only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(model_a)))
)
G_a = build_graph(data_a; time_step = 1)
add_hvdc_edges!(G_a, sys_a, res_vars_a, data_a; time_step = 1)
branches_a = collect(PSY.get_components(PSY.ACBranch, sys_a))
C_a = find_circular_flows(G_a, data_a, branches_a)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: buses=$(c.bus_numbers), branches=$(c.branches), min_flow=$(minimum(c.branch_flows)) MW")
end

println("=== Scenario B: without circular flow (negative RE cost_sign) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0)
model_b = build_rts_model_with_quadratic_losses(sys_b)
PSI.solve!(model_b)
results_b = PSI.OptimizationProblemResults(model_b)
res_vars_b = results_b.variable_values
data_b = PSI.get_power_flow_data(
    only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(model_b)))
)
G_b = build_graph(data_b; time_step = 1)
add_hvdc_edges!(G_b, sys_b, res_vars_b, data_b; time_step = 1)
branches_b = collect(PSY.get_components(PSY.ACBranch, sys_b))
C_b = find_circular_flows(G_b, data_b, branches_b)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), branches=$(c.branches), min_flow=$(minimum(c.branch_flows)) MW")
end
