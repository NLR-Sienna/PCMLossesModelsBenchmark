# End-to-end build validation of the CircularFlows module with quadratic losses on the RTS
# test system.  Two scenarios demonstrate the same setup as rts_example.jl but with the
# quadratic loss term  loss = -Σₖ Rₖ(Σⱼ PTDFₖⱼ Pⱼ)²  added to the copper-plate balance.
#
# Scenario A: renewable cost > 0  → HVDC loop is profitable → circular flow expected
# Scenario B: HVDC disabled       → no inter-area loop possible → no circular flow
#
# NOTE: Solving the nonconvex MIQP requires Gurobi with NonConvex=2.  This script only
# validates that build! succeeds with the quadratic loss constraints in place.

using Pkg
this_path = @__DIR__
Pkg.activate(joinpath(this_path, "..", ".."))

using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using HiGHS
using Graphs
using SimpleWeightedGraphs
using Dates
using DataFrames

include(joinpath(this_path, "mapped_indices.jl"))
include(joinpath(this_path, "circular_flows.jl"))
include(joinpath(this_path, "../../Systems/RTS/build_rts.jl"))
include(joinpath(this_path, "../build_models.jl"))
include(joinpath(this_path, "../utils.jl"))

# ---- Scenario A: negative renewable costs → circular flow expected when solved ----
println("=== Scenario A: with circular flow (negative RE costs) ===")
sys_a = build_rts_system()
set_renewable_costs!(sys_a, 1.0) # Positive number is added negative to the objective function
model_a = build_rts_model_with_quadratic_losses(sys_a)
println("  Model A built successfully with quadratic losses ✓")
# To detect circular flows after solving with Gurobi (NonConvex=2):
#   PSI.solve!(model_a)
#   results_a = PSI.OptimizationProblemResults(model_a)
#   res_vars_a = PSI.read_variables(results_a)
#   data_a = PSI.get_power_flow_data(
#       only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(model_a))))
#   G_a = build_graph(data_a; time_step=1)
#   add_hvdc_edges!(G_a, sys_a, res_vars_a, data_a; time_step=1)
#   branches_a = collect(PSY.get_components(PSY.ACBranch, sys_a))
#   C_a = find_circular_flows(G_a, data_a, branches_a)
#   @assert length(C_a) > 0 "Scenario A FAILED: expected ≥1 circular flow, got 0"

# ---- Scenario B: positive renewable costs → circular flow not expected when solved ----
println("=== Scenario B: without circular flow (HVDC disabled) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0) # Negative number is added positive to the objective function, but HVDC loop incentive is removed by setting HVDC losses to 0 and allowing free reactive power
model_b = build_rts_model_with_quadratic_losses(sys_b)
println("  Model B built successfully with quadratic losses ✓")
# To detect circular flows after solving with Gurobi (NonConvex=2):
#   PSI.solve!(model_b)
#   results_b = PSI.OptimizationProblemResults(model_b)
#   res_vars_b = PSI.read_variables(results_b)
#   data_b = PSI.get_power_flow_data(
#       only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(model_b))))
#   G_b = build_graph(data_b; time_step=1)
#   add_hvdc_edges!(G_b, sys_b, res_vars_b, data_b; time_step=1)
#   branches_b = collect(PSY.get_components(PSY.ACBranch, sys_b))
#   C_b = find_circular_flows(G_b, data_b, branches_b)
#   @assert length(C_b) == 0 "Scenario B FAILED: expected 0 circular flows, got $(length(C_b))"

println("\nBoth models built with quadratic losses.")
println("Solve with Gurobi (NonConvex=2) to run circular flow detection.")
