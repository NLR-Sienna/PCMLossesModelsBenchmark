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
# =============================================================================
# PART 5 – FLOW CANCELLING WITH QUADRATIC LOSS APPROXIMATION
# =============================================================================
# Combines the flow-cancelling investment model (Part 4) with the quadratic
# loss formulation (Part 3).  Requires an NLP-capable solver because of the
# quadratic loss constraints.

model_fc_quad = build_model_with_flow_canceling_and_quadratic_losses(
    sys_inv;
    optimizer = optimizer_with_attributes(Gurobi.Optimizer),
)

container_quad = model_fc_quad.internal.container

# #Relax binary for IPOPT to solve
# opt_var = all_variables(container_quad.JuMPmodel)
# for v in opt_var
#     if is_binary(v)
#         unset_binary(v); set_upper_bound(v, 1)
#         set_lower_bound(v, 0);      
#     end
# end
# ──► build_model_with_flow_canceling_and_quadratic_losses
#   [FlowCancelling/build_models.jl:666-714]
#
# Steps 1-8 are identical to build_model_with_flow_canceling_terms above.
#
# Then three additional steps add the quadratic loss approximation:
#
# Step 9: _fc_add_loss_variables!(model)
#   [FlowCancelling/build_models.jl:575-595]
#   Identical role to add_current_loss_variables! in the losses path:
#   adds LineLossTotalApproximation[ref_bus, t]  continuous variable.
#
# Step 10: _fc_add_loss_to_copperplate_balance!(model, loss_var)
#   [FlowCancelling/build_models.jl:600-609]
#   Modifies CopperPlateBalanceConstraint: Σ injection + loss_var = 0
#   (quadratic variant: RHS stays 0, loss value set by Step 11)
#
# Step 11: _fc_add_quadratic_loss_constraints!(model, sys, ptdf)
#   [FlowCancelling/build_models.jl:616-651]
#   Adds LineLossConstraintApproximation (quadratic, no voltage scaling):
#   loss_var[t] = -Σ_k R[k] * FC_flow[k,t]²
#   where FC_flow[k,t] is the PTDFBranchFlowWithFC expression.

# PTDFBranchFlowWithFC expressions are added to the container during the build
# step above (before solve!).  Explore both PTDFBranchFlow and
# PTDFBranchFlowWithFC directly via model_fc_quad.internal.container.expressions.
#
# The expressions dict maps ExpressionKey → DenseAxisArray{JuMP.AffExpr}.
# Printing a JuMP.AffExpr entry shows the symbolic structure: PTDF coefficients
# multiplied by injection VariableRefs (PTDFBranchFlow), plus BranchCancellingFlow
# VariableRefs with shift-factor coefficients (PTDFBranchFlowWithFC).

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

# Read the optimal FC-corrected branch flows.
# Because should_write_resulting_value(PTDFBranchFlowWithFC) = true, the
# solved expression values are stored in the result alongside variables.
# PTDFBranchFlowWithFC[branch, t] incorporates the investment-dependent
# flow-cancelling correction, so the values reflect the as-built topology.
# fc_flows_quad = read_expression(res_fc_quad, "PTDFBranchFlowWithFC__Line")
# println("FC-corrected branch flows at optimality (Line):\n", fc_flows_quad)

#ychen
println("Solved Losses,",
value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))
#Losses are negatve: Gen+Loss=Load
#
#Calculate losses using r*I^2
fc_line = PSI.get_expression(container_quad, PTDFBranchFlowWithFC(), PSY.Line)
line_names = axes(fc_line, 1)
R_line = Dict(get_name(l) => get_r(l) for l in get_components(get_available, PSY.Line, sys_inv))

other_branch_type_data = []
for T in (PSY.PhaseShiftingTransformer, PSY.TapTransformer, PSY.Transformer2W)
    key = InfrastructureSystems.Optimization.ExpressionKey{PTDFBranchFlowWithFC, T}("")
    if !haskey(container_quad.expressions, key)
        continue
    end
    fc = container_quad.expressions[key]
    R_dict = Dict(get_name(b) => get_r(b) for b in get_components(get_available, T, sys_inv))
    push!(other_branch_type_data, (fc, R_dict))
end

cal_loss=Dict()
for t in axes(fc_line, 2)
    cal_loss[t] = 
        -sum(R_line[name] * value(fc_line[name, t])^2 for name in line_names) -
        sum(
            R_dict[name] * value(fc[name, t])^2
            for (fc, R_dict) in other_branch_type_data
            for name in axes(fc, 1);
            init = 0.0,
        )
end        
#julia> cal_loss 
#Dict{Any, Any} with 2 entries:
#  2 => -0.0585694
#  1 => -0.0494142

#=
uc1=model_fc_quad.internal.container
open("model_fc_quad.txt","w") do io
    redirect_stdout(io) do
        println(objective_function(uc1.JuMPmodel))
        for k in all_constraints(uc1.JuMPmodel,; include_variable_in_set_constraints = true)           
            println(name(k),",",k) 
        end    
    end
end
=#
#set candidate_line_2 investment to 1 and candidate_line_1 to 0 base on the relaxed solution
#value.(container_quad.variables[PSI.VariableKey{BranchInvestmentVariable, Line}("")])
#1-dimensional DenseAxisArray{Float64,1,...} with index sets:
#    Dimension 1, ["candidate_line_2", "candidate_line_1"]
#And data, a 2-element Vector{Float64}:
# 0.9872490778825471
# 1.9270021486612987e-8
#julia> value.(container_quad.variables[PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}("")])
#1-dimensional DenseAxisArray{Float64,1,...} with index sets:
#    Dimension 1, ["candidate_thermal_1", "candidate_thermal_2"]
#And data, a 2-element Vector{Float64}:
#  0.9999995441878543
# -6.641140368302735e-9
#fix(container_quad.variables[PSI.VariableKey{BranchInvestmentVariable, Line}("")]["candidate_datacenter_line_1"],0,force=true)
fix(container_quad.variables[PSI.VariableKey{BranchInvestmentVariable, Line}("")]["candidate_datacenter_line_2_345kV"],1,force=true)
fix(container_quad.variables[PSI.VariableKey{BranchInvestmentVariable, Line}("")]["candidate_datacenter_line_1_500kV"],0,force=true)
#fix(container_quad.variables[PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}("")]["candidate_thermal_1"],1,force=true)
fix(container_quad.variables[PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}("")]["candidate_thermal_2"],1,force=true)

solve!(model_fc_quad)
value.(container_quad.variables[PSI.VariableKey{BranchInvestmentVariable, Line}("")])
value.(container_quad.variables[PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}("")])

println("Solved Losses,",
value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))

#julia> println("Solved Losses,",
#       value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))
#Solved Losses,2-dimensional DenseAxisArray{Float64,2,...} with index sets:
#    Dimension 1, [4]
#    Dimension 2, 1:2
#And data, a 1×2 Matrix{Float64}:
# -0.047405713992274  -0.05802386784571255



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

#add_candidate_line_data_without_parallel!(sys1)
add_datacenter_data!(sys1)
#add_candidate_datacenter_line_data!(sys1, sc)
add_candidate_datacenter_line_data_kv!(sys1, 230, 345)

unbuilt_line=get_component(Line,sys1,"candidate_datacenter_line_2_345kV")
#unbuilt_gen=get_component(Generator,sys1,"candidate_thermal_2")

set_available!(unbuilt_line,false)
#set_available!(unbuilt_gen,false)

#transform_single_time_series!(sys1, Hour(2), Hour(2))

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
    optimizer = Gurobi.Optimizer,
    name = "UC_QuadLoss",
    store_variable_names = true,
)

build!(model_bench; output_dir = mktempdir())
    # --- Quadratic loss approximation ---
loss_var = add_current_loss_variables!(model_bench)

    # Step 2: Modify nodal balance to account for losses (without fixed RHS)
add_quadratic_current_loss_to_copperplate_balance!(model_bench, loss_var)

    # Step 3: Add quadratic loss approximation constraints (P = I²R formulation)
add_current_loss_constraint_quadratic_approximation_no_voltage!(model_bench, sys1, ptdf)
solve!(model_bench)

container_quad = model_bench.internal.container

println("Solved Losses,",
value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))

#=
julia> println("Solved Losses,",                                                                                                
       value.(container_quad.variables[InfrastructureSystems.Optimization.VariableKey{LineLossTotalApproximation, System}("")]))
Solved Losses,2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, [4]
    Dimension 2, 1:2
And data, a 1×2 Matrix{Float64}:
 -0.04740571263845791  -0.05802386747719086
=#

## Flow cancelling lossless model:
# sys_inv = build_matpower_5bus_with_updated_lines()
# transform_single_time_series!(sys_inv, Hour(2), Hour(2))
# # Comment an additional phase shifting transformer to avoid an issue with different parallel types in PSI
# set_available!(get_component(PhaseShiftingTransformer, sys_inv, "bus-3-bus-4-i_5"), false)
# sc = "parameters" #scenario = "parameters" # "parameters", "costs", or "both" (see build_5bus_datacenter_update.jl for details)

# # Add candidate thermal generators (flagged with ext["is_candidate"] = true)
# # ──► candidate_projects_data  [Systems/5bus/build_5bus.jl:44-99]
# candidate_gens = candidate_projects_data(sys_inv)
# for gen in candidate_gens
#     add_component!(sys_inv, gen)
# end

# # Add candidate transmission lines using the "no-parallel" topology:
# # each existing line is split into two segments with an intermediate bus so
# # that the candidate line shares the arc but not the physical conductor.
# # ──► add_candidate_line_data_without_parallel!  [Systems/5bus/build_5bus.jl:268-277]
# #add_candidate_line_data_without_parallel!(sys_inv)
# add_datacenter_data!(sys_inv)
# add_candidate_datacenter_line_data_kv!(sys_inv, 230, 345)
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