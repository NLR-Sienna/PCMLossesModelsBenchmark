using Pkg
Pkg.activate(".")

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
using Xpress
using Logging

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

include("Systems/5bus/build_5bus.jl")                      # 5-bus test-case helpers
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
add_candidate_line_data_without_parallel!(sys_inv)

# Now build the model with flow-cancelling terms.
model_fc = build_model_with_flow_canceling_terms(sys_inv)
# ──► build_model_with_flow_canceling_terms  [FlowCancelling/build_models.jl:497-560]
#
# Step 1: build base PTDF DecisionModel (identical to Part 1 structure)
#
# Step 2: add_branch_investment_variables!(model, candidate_lines, Line)
#   [FlowCancelling/build_models.jl:158-178]
#   Adds binary variable z_k ∈ {0,1} per candidate line
#   (PSI.add_variable_container! → BranchInvestmentVariable)
#
# Step 3: add_branch_cancelling_flow_variables!(model, candidate_lines, Line)
#   [FlowCancelling/build_models.jl:187-208]
#   Adds continuous variable v_k[t] per (candidate line, time step)
#   (PSI.add_variable_container! → BranchCancellingFlowVariable)
#
# Step 4: add_bigM_linking_constraints!(model, candidate_lines, Line, z_var, v_var)
#   [FlowCancelling/build_models.jl:222-256]
#   Adds BigMConstraint (ub and lb):
#     v_k[t] ≤  M_max * (1 - z_k)
#     v_k[t] ≥ -M_max * (1 - z_k)
#   When z_k = 1 (line built): v_k[t] is forced to 0  → no cancellation needed
#   When z_k = 0 (line not built): v_k[t] free in [-M, M] → cancels spurious flow
#
# Step 5: add_shift_terms_to_existing_line_constraints! (called for Line and
#         PhaseShiftingTransformer)  [FlowCancelling/build_models.jl:266-298]
#   For each existing branch l and each candidate line k, appends
#   Δ_{l,k} * v_k[t] to the FlowRateConstraint ub and lb:
#     set_normalized_coefficient(FlowRateConstraint_ub[l,t], v_var[k,t], shift)
#     set_normalized_coefficient(FlowRateConstraint_lb[l,t], v_var[k,t], shift)
#   where Δ_{l,k} = PTDF[from_k, l] - PTDF[to_k, l]  (shift factor, eq. 3)
#
# Step 6: add_shift_terms_to_candidate_line_constraints!
#   [FlowCancelling/build_models.jl:312-369]
#   Modifies FlowRateConstraint bounds for candidate lines (eq. 18):
#     UB: flow[k,t] - rating_k * z_k + (Δ_{k,k} - 1) * v_k[t]
#         + Σ_{j≠k} Δ_{k,j} * v_j[t]  ≤  0
#     LB: -flow[k,t] - rating_k * z_k - (Δ_{k,k} - 1) * v_k[t]
#         - Σ_{j≠k} Δ_{k,j} * v_j[t]  ≤  0
#   When z_k = 0: flow[k,t] = 0  (line not in service)
#   When z_k = 1: flow[k,t] ∈ [-rating_k, +rating_k]  (normal operation)
#
# Step 7: add_candidate_line_investment_costs!(model, z_var)
#   [FlowCancelling/build_models.jl:387-405]
#   Reads project_cost from ext dict of each candidate line;
#   appends  cost_k * z_k  to the JuMP objective via
#   JuMP.add_to_expression!(objective_function(jump_model), cost, z_var[name])
#
# Step 8: add_candidate_generation_investment_constraints!(model, ThermalStandard)
#   [FlowCancelling/build_models.jl:420-480]
#   For each generator with ext["is_candidate"] = true:
#     - Adds binary variable x_g  (GenerationInvestmentVariable)
#     - Adds constraint: p_g[t] ≤ p_max_g * x_g  (GenerationInvestmentConstraint)
#     - Appends project_cost_g * x_g to objective

solve!(model_fc)

res_fc = OptimizationProblemResults(model_fc)
println("=== Flow-cancelling model (lossless) solved ===")
println("Objective: ", JuMP.objective_value(model_fc.internal.container.JuMPmodel))

# Read investment decisions
inv_lines = read_variable(res_fc, PSI.VariableKey{BranchInvestmentVariable, Line}(""))
inv_gens  = read_variable(res_fc, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions:\n", inv_lines)
println("Generation investment decisions:\n", inv_gens)

# Explore both PTDFBranchFlow and PTDFBranchFlowWithFC expressions directly
# from the JuMP model container.  Both are registered during build_model_with_-
# flow_canceling_terms — PTDFBranchFlowWithFC is available even in the lossless
# model so the reader can inspect the FC correction without running the NLP.
#
# model_fc.internal.container.expressions is a Dict{ExpressionKey, DenseAxisArray}
# where each array entry is a JuMP.AffExpr — a symbolic linear combination of
# VariableRefs that can be printed before any solve to see which variables
# appear and with what coefficients.
ptdf_key_fc = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlow,       Line}("")
fc_key_fc   = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlowWithFC, Line}("")

ptdf_exprs_fc = model_fc.internal.container.expressions[ptdf_key_fc]
fc_exprs_fc   = model_fc.internal.container.expressions[fc_key_fc]
first_branch  = first(axes(ptdf_exprs_fc, 1))

println("PTDFBranchFlow[\"$first_branch\", 1] as JuMP.AffExpr (PTDF·injection only):")
println("  ", ptdf_exprs_fc[first_branch, 1])
println("PTDFBranchFlowWithFC[\"$first_branch\", 1] as JuMP.AffExpr (+ FC correction terms):")
println("  ", fc_exprs_fc[first_branch, 1])

# Read the same flows at the optimal dispatch as numerical values.
# PTDFBranchFlow[branch, t] = Σⱼ PTDF[branch,j] · injection[j,t]
ptdf_flows_line = read_expression(res_fc, "PTDFBranchFlow__Line")
println("PTDF branch flows at optimality (Line):\n", ptdf_flows_line)

#Ychen error PTDFBranchFlowWithFC__Line is not stored
# can be retrieved from the optimization 
#fc_flows_line = read_expression(res_fc, "PTDFBranchFlowWithFC__Line")
#println("FC-corrected branch flows at optimality (Line):\n", fc_flows_line)

#flow without flow cancelling terms
println(value.(ptdf_exprs_fc))

#flow with flow cancelling terms
println(value.(fc_exprs_fc))
#Ychen This shows candidate_line_1 flow =0, matching candidate_line_1 inv_lines inv_line=0
#=
julia> println(value.(fc_exprs_fc))
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, ["bus-1-bus-10-i_3", "bus-1-bus-2-i_1", "bus-1-bus-4-i_2_segment_1", "bus-1-bus-4-i_2_segment_2", "bus-2-bus-3-i_4", "bus-4-bus-10-i_7_segment_1", "bus-4-bus-10-i_7_segment_2", "candidate_line_1", "candidate_line_2"]
    Dimension 2, 1:2
And data, a 9×2 Matrix{Float64}:
 -2.7539364126035197     -4.001055883691279
  2.505593680282442       2.7010326362317363
  1.3741713661605393      1.1105883247297712
  1.3741713661605393      1.1105883247297712
 -0.4944063197175584      0.5314279142317362
 -1.9999999999999996     -1.998944116308719
 -1.9999999999999996     -1.998944116308719
  4.440892098500626e-16  -2.220446049250313e-16
  1.3741713661605393      1.1105883247297712
=#  

#####################################################################
###Ychen validate through setting candidate_line_1 unavailable
sys1 = build_matpower_5bus_with_updated_lines()
transform_single_time_series!(sys1, Hour(2), Hour(2))
# Comment an additional phase shifting transformer to avoid an issue with different parallel types in PSI
set_available!(get_component(PhaseShiftingTransformer, sys1, "bus-3-bus-4-i_5"), false)

# Add candidate thermal generators (flagged with ext["is_candidate"] = true)
# ──► candidate_projects_data  [Systems/5bus/build_5bus.jl:44-99]
candidate_gens = candidate_projects_data(sys1)
for gen in candidate_gens
    add_component!(sys1, gen)
end

add_candidate_line_data_without_parallel!(sys1)

unbuilt_line=get_component(Line,sys1,"candidate_line_1")
unbuilt_gen=get_component(Generator,sys1,"candidate_thermal_2")

set_available!(unbuilt_line,false)
set_available!(unbuilt_gen,false)

transform_single_time_series!(sys1, Hour(2), Hour(2))

ptdf = PTDF(sys1)

network_model = NetworkModel(PTDFPowerModel; PTDF_matrix = ptdf, use_slacks = true)

template = ProblemTemplate(network_model)
set_device_model!(template, ThermalStandard, ThermalDispatchNoMin)
set_device_model!(template, Line, StaticBranch)
set_device_model!(template, PhaseShiftingTransformer, StaticBranch)
set_device_model!(template, PowerLoad, StaticPowerLoad)

model_bench = DecisionModel(
    template,
    sys1;
    optimizer = DEFAULT_MILP_OPTIMIZER, #Xpress.Optimizer,
    name = "UC",
    store_variable_names=true,
)

build!(model_bench; output_dir = mktempdir())

solve!(model_bench)

ptdf_exprs_bench = model_bench.internal.container.expressions[ptdf_key_fc]

println(value.(ptdf_exprs_bench))
#ychen The values match value.(fc_exprs_fc)
#=
julia> println(value.(ptdf_exprs_bench))
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, ["bus-1-bus-10-i_3", "bus-1-bus-2-i_1", "bus-1-bus-4-i_2_segment_1", "bus-1-bus-4-i_2_segment_2", "bus-2-bus-3-i_4", "bus-4-bus-10-i_7_segment_1", "bus-4-bus-10-i_7_segment_2", "candidate_line_2"]
    Dimension 2, 1:2
And data, a 8×2 Matrix{Float64}:
 -2.753936412603519   -4.001055883691281
  2.5055936802824403   2.7010326362317367
  1.3741713661605393   1.110588324729772
  1.3741713661605393   1.110588324729772
 -0.4944063197175601   0.5314279142317366
 -2.0                 -1.9989441163087203
 -2.0                 -1.9989441163087203
  1.3741713661605393   1.110588324729772
  =#

  #================