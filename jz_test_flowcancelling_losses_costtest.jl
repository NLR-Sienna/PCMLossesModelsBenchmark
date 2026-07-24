# using Pkg
# Pkg.activate(".")

using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using LinearAlgebra
using HiGHS
using Ipopt
using Dates
#using Xpress
using Logging
using Gurobi
import PowerSystems as PSY
import PowerSimulations as PSI
import PowerSystemCaseBuilder as PSB

# Project-local scripts.  Order matters: build_models.jl defines constants and
# structs that run_models.jl and utils.jl depend on.
include("SiennaScripts/build_models.jl")   # loss approximation builders
include("SiennaScripts/utils.jl")          # post-processing helpers (loss factors, FND, …)
include("SiennaScripts/run_models.jl")     # iterative solve loops

# =============================================================================
# PART 4 – FLOW CANCELLING (transmission expansion, lossless)
# =============================================================================
include("Systems/5bus/ac_line_expansion_model_example.jl")
include("Systems/5bus/build_5bus_datacenter_update.jl")                      # 5-bus test-case helpers
include("SiennaScripts/FlowCancelling/build_models.jl")    # flow-cancelling builders

# Flow cancelling enables an investment model where candidate lines can be
# "connected" (binary z_k = 1) or ignored (z_k = 0).  When z_k = 0 the
# candidate line still appears in the PTDF matrix (it was there at build!
# time), so its flow would normally affect other branches.  The flow-cancelling
# variables v_k[t] ≈ z_k * flow[k,t] are used to subtract out those spurious
# flows when the line is not built.

# Build the 5-bus investment system: existing generators + candidate generators
# + existing lines + candidate lines (parallel topology via intermediate buses)
sys_inv = build_matpower_5bus_with_updated_lines()
transform_single_time_series!(sys_inv, Hour(2), Hour(2))
# Comment an additional phase shifting transformer to avoid an issue with different parallel types in PSI
set_available!(get_component(PhaseShiftingTransformer, sys_inv, "bus-3-bus-4-i_5"), false)
sc = "both" #scenario = "parameters" # "parameters", "costs", or "both" (see build_5bus_datacenter_update.jl for details)

# Add candidate thermal generators (flagged with ext["is_candidate"] = true)
# ──► candidate_projects_data  [Systems/5bus/build_5bus.jl:44-99]
candidate_gens = candidate_projects_data(sys_inv)
for gen in candidate_gens
    add_component!(sys_inv, gen)
end

# Add candidate transmission lines using the "no-parallel" topology:
# each existing line is split into two segments with an intermediate bus so
# that the candidate line shares the arc but not the physical conductor.
# ──► add_candidate_line_data_without_parallel!  [Systems/5bus/build_5bus.jl:268-277]
#add_candidate_line_data_without_parallel!(sys_inv)
add_datacenter_data!(sys_inv)
#add_candidate_datacenter_line_data!(sys_inv, sc)
add_candidate_datacenter_line_data_kv!(sys_inv, 230, 345)
#add_candidate_datacenter_line_data_cost!(sys_inv, 0)

model_fc_quad = build_model_with_flow_canceling_and_quadratic_losses(
    sys_inv;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
)

container_quad = model_fc_quad.internal.container

ptdf_key_quad = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlow,       Line}("")
fc_key_quad   = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlowWithFC, Line}("")

println("Expression keys registered in the container:")
for k in keys(container_quad.expressions); println("  ", k); end

ptdf_exprs_quad = container_quad.expressions[ptdf_key_quad]
fc_exprs_quad   = container_quad.expressions[fc_key_quad]
first_line      = first(axes(fc_exprs_quad, 1))

println("PTDFBranchFlow[\"$first_line\", 1] (PTDF·injection, no FC correction):")
println("  ", ptdf_exprs_quad[first_line, 1])
println("PTDFBranchFlowWithFC[\"$first_line\", 1] (adds BranchCancellingFlow correction terms):")
println("  ", fc_exprs_quad[first_line, 1])

solve!(model_fc_quad)

# This will fail since we are using Ipopt that does not support binary variables.
# However, Gurobi can be used to solve MINLP problems.

res_fc_quad = OptimizationProblemResults(model_fc_quad)
println("=== Flow-cancelling + quadratic losses model solved ===")
println("Objective: ", JuMP.objective_value(model_fc_quad.internal.container.JuMPmodel))

inv_lines_quad = read_variable(res_fc_quad, PSI.VariableKey{BranchInvestmentVariable, Line}(""))
inv_gens_quad  = read_variable(res_fc_quad, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions (with losses):\n", inv_lines_quad)
println("Generation investment decisions (with losses):\n", inv_gens_quad)

println("Solved Losses,",
value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))
#GRBterminate(backend(model_fc_quad))
#GRBfreemodel(backend(model_fc_quad))

model_fc_lossless = build_model_with_flow_canceling_terms(
    sys_inv;
)

container_lossless = model_fc_lossless.internal.container

ptdf_key_lossless = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlow,       Line}("")
fc_key_lossless   = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlowWithFC, Line}("")

println("Expression keys registered in the container:")
for k in keys(container_lossless.expressions); println("  ", k); end

ptdf_exprs_lossless = container_lossless.expressions[ptdf_key_lossless]
fc_exprs_lossless   = container_lossless.expressions[fc_key_lossless]
first_line      = first(axes(fc_exprs_lossless, 1))

println("PTDFBranchFlow[\"$first_line\", 1] (PTDF·injection, no FC correction):")
println("  ", ptdf_exprs_lossless[first_line, 1])
println("PTDFBranchFlowWithFC[\"$first_line\", 1] (adds BranchCancellingFlow correction terms):")
println("  ", fc_exprs_lossless[first_line, 1])

solve!(model_fc_lossless)

# This will fail since we are using Ipopt that does not support binary variables.
# However, Gurobi can be used to solve MINLP problems.

res_fc_lossless = OptimizationProblemResults(model_fc_lossless)
println("=== Flow-cancelling model solved ===")
println("Objective: ", JuMP.objective_value(model_fc_lossless.internal.container.JuMPmodel))

inv_lines_lossless = read_variable(res_fc_lossless, PSI.VariableKey{BranchInvestmentVariable, Line}(""))
inv_gens_lossless  = read_variable(res_fc_lossless, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions (with losses):\n", inv_lines_lossless)
println("Generation investment decisions (with losses):\n", inv_gens_lossless)

# Read the optimal FC-corrected branch flows.
# Because should_write_resulting_value(PTDFBranchFlowWithFC) = true, the
# solved expression values are stored in the result alongside variables.
# PTDFBranchFlowWithFC[branch, t] incorporates the investment-dependent
# flow-cancelling correction, so the values reflect the as-built topology.
# fc_flows_quad = read_expression(res_fc_quad, "PTDFBranchFlowWithFC__Line")
# println("FC-corrected branch flows at optimality (Line):\n", fc_flows_quad)

#ychen
# println("Solved Losses,",
# value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))