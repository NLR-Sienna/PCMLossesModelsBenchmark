# =============================================================================
# TUTORIAL: Transmission Loss Approximations and Flow Cancelling
# =============================================================================
#
# PURPOSE
# -------
# This script is a guided tour of the two families of methods added on top of
# PowerSimulations.jl (PSI) in this project:
#
#   1. LOSS APPROXIMATIONS – incorporating transmission losses into a standard
#      PTDF-based unit commitment / economic dispatch problem.
#      Files: SiennaScripts/build_models.jl, SiennaScripts/run_models.jl,
#             SiennaScripts/utils.jl
#
#   2. FLOW CANCELLING – a topology-expansion investment model where candidate
#      lines are selectively connected / disconnected via binary variables and
#      big-M constraints (Caramanis et al., IEEE PES 2016).
#      File: SiennaScripts/FlowCancelling/build_models.jl
#
# KEY MENTAL MODEL
# ----------------
# PSI builds a JuMP model from a ProblemTemplate and a PowerSystems.jl System.
# The methods in this project do NOT change that build step.  Instead, they
# reach into the already-built `DecisionModel` and surgically modify its JuMP
# model by:
#   * adding new JuMP variables via  PSI.add_variable_container!
#   * adding new JuMP constraints via PSI.add_constraints_container! + @constraint
#   * modifying existing constraint coefficients / RHS via
#     set_normalized_coefficient / set_normalized_rhs
#   * modifying the objective via JuMP.add_to_expression!
#
# Comments throughout this script point to the exact functions that implement
# each of these modifications, so you can jump straight to the source.
#
# HOW TO USE THIS SCRIPT
# ----------------------
# Run it top-to-bottom once to see the full workflow.
# After that, use the function-name annotations as a map into the source:
# each comment of the form  ──► FunctionName  [file:line-range]
# tells you exactly where the interesting code lives.

# =============================================================================
# 0. ENVIRONMENT SETUP
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using HiGHS
using Ipopt
using Dates
using Logging

import PowerSystems as PSY
import PowerSimulations as PSI
import PowerSystemCaseBuilder as PSB

# Project-local scripts.  Order matters: build_models.jl defines constants and
# structs that run_models.jl and utils.jl depend on.
include("SiennaScripts/build_models.jl")   # loss approximation builders
include("SiennaScripts/utils.jl")          # post-processing helpers (loss factors, FND, …)
include("SiennaScripts/run_models.jl")     # iterative solve loops

include("Systems/5bus/build_5bus.jl")                      # 5-bus test-case helpers
include("SiennaScripts/FlowCancelling/build_models.jl")    # flow-cancelling builders

# =============================================================================
# PART 1 – LOSSLESS BASELINE (standard PSI PTDF model)
# =============================================================================
# We start with the standard PowerSimulations PTDF-based unit commitment.
# No losses, no new variables – just the textbook formulation as PSI builds it.

sys_uc = PSB.build_system(PSITestSystems, "c_sys5_uc";
    skip_serialization = true,
    add_single_time_series = true,
)
set_available!(get_component(TwoTerminalGenericHVDCLine, sys_uc, "DC1"), false)
transform_single_time_series!(sys_uc, Hour(1), Hour(1))

# ──► build_ptdf_model_without_losses  [SiennaScripts/build_models.jl:201-217]
#
# Internally calls make_ptdf_model_without_losses (lines 69-123), which:
#   1. Builds a ProblemTemplate with NetworkModel(PTDFPowerModel):
#        - PTDF_matrix: pre-computed PTDF object
#        - use_slacks:  slack variables on every constraint for infeasibility diagnosis
#        - duals:       requests dual values for CopperPlateBalanceConstraint
#        - power_flow_evaluation: runs an AC power flow AFTER each solve to compute
#          loss factors (PowerFlowLossFactors, PowerFlowVoltageMagnitude, etc.)
#   2. Adds device formulations via set_device_model! for each entry in device_models
#   3. Creates the DecisionModel and calls build!(model; output_dir = mktempdir())
#
# After build!, the JuMP model contains:
#   Copper-plate balance:   Σ_i injection[i,t]  = 0   ∀ t
#   Branch flow constraints: flow[k,t] ∈ [-rating_k, +rating_k]  ∀ k, t
#   where flow[k,t] = Σ_i PTDF[k,i] * injection[i,t]
model_lossless = build_ptdf_model_without_losses(sys_uc)

solve!(model_lossless)
res_lossless = OptimizationProblemResults(model_lossless)

println("=== Lossless model solved ===")
println("Objective: ", res_lossless.optimizer_stats[1, "objective_value"])

# At this point res_lossless carries auxiliary variables from the AC power flow:
#   read_aux_variable(res_lossless, "PowerFlowLossFactors__ACBus")       – ∂Loss/∂P per bus
#   read_aux_variable(res_lossless, "PowerFlowVoltageMagnitude__ACBus")  – |V| per bus
#   read_aux_variable(res_lossless, "PowerFlowBranchActivePowerFromTo__Line")  – branch MW
#
# These are the raw inputs for all loss-approximation methods below.

# =============================================================================
# PART 2 – ITERATIVE LINEAR LOSS APPROXIMATION
# =============================================================================
# Real transmission networks have losses.  The standard PTDF model ignores them
# because they introduce a nonlinearity (P_loss ∝ I² ∝ flow²).
#
# The linear approach replaces the nonlinear loss with a first-order Taylor
# expansion around the previous operating point, then re-solves.  Repeating
# this until the dispatch stops changing gives a converged estimate.

# ──► run_iterative_linear_loss_model  [SiennaScripts/run_models.jl:170-226]
#
# Algorithm:
#   iter 0:  solve lossless model   (run_lossless_model, lines 25-59)
#   iter k:  solve linear-loss model (run_linear_loss_model, lines 61-124)
#            stop when max|P_gen_new - P_gen_old| < error_tol
#
# Each call to run_linear_loss_model calls build_ptdf_model_with_linear_losses
# (lines 809-838 in build_models.jl), which does four things after build!:
#
#   STEP A – get_bus_loss_factors(res_old)         [utils.jl:73-89]
#     Reads PowerFlowLossFactors__ACBus from the previous solve.
#     Returns matrix loss_factors[bus, t] = ∂TotalLoss/∂Injection[bus,t]
#     (computed by PSI's AC power flow evaluation after each solve)
#
#   STEP B – get_total_AC_loss(res_old)            [utils.jl:418-473]
#     Sums P_from_to + P_to_from across all branches to get the "true" AC
#     losses in MW per time step.  This is the RHS offset for the balance eq.
#
#   STEP C – update_copperplate_loss_approximation!(model, loss_factors,
#                                                   total_loss_est, injection_old)
#     [build_models.jl:620-634]  Three sub-steps:
#
#     C1. add_current_loss_variables!              [build_models.jl:234-260]
#           PSI.add_variable_container! → JuMP.@variable  per (ref_bus, t)
#           New variable: LineLossTotalApproximation[ref_bus, t]  (continuous, unbounded)
#
#     C2. add_current_loss_to_copperplate_balance! [build_models.jl:286-311]
#           set_normalized_coefficient(CopperPlateBalanceConstraint[t], loss_var[t], 1)
#           set_normalized_rhs(CopperPlateBalanceConstraint[t],
#                              rhs + total_loss_est[t] / base_power)
#           Before: Σ injection[i,t] = 0
#           After:  Σ injection[i,t] + loss_var[t] = total_loss_est[t] / base_power
#           → forces the optimizer to procure enough generation to cover estimated losses
#
#     C3. add_current_loss_constraint_approximation! [build_models.jl:383-421]
#           Adds a new constraint (LineLossConstraintApproximation):
#           loss_var[t] = Σ_i (injection_old[i,t] - injection[i,t]) * loss_factors[i,t]
#           This is the Taylor expansion:
#           ΔLoss ≈ Σ_i (∂Loss/∂P_i) * ΔP_i
#
#   STEP D – update_transmission_constraints_with_losses!(model, res_old, sys, ptdf)
#     [build_models.jl:722-772]
#     Adjusts the RHS of FlowRateConstraint (branch flow limits) to account for
#     loss-induced loading via the Fictitious Nodal Demand (FND) approach:
#       - get_fictitious_nodal_demand_by_loss(res_old, sys)  [utils.jl:806-881]
#         splits each branch's loss equally between its two endpoint buses
#       - For each branch k: new_RHS = old_RHS + Σ_i PTDF[k,i] * FND[i,t]

res_linear = run_iterative_linear_loss_model(sys_uc; max_iter = 10, error_tol = 1e-3)

println("=== Iterative linear loss model converged ===")
println("Objective: ", res_linear.optimizer_stats[1, "objective_value"])

# =============================================================================
# PART 3 – QUADRATIC LOSS APPROXIMATION (single-shot, NLP)
# =============================================================================
# Instead of linearizing, model losses as P_loss = Σ_k R_k * (PTDF_k * P)²
# This is physically accurate (Ohm's law) but produces a quadratic program
# requiring an NLP solver (Ipopt).  Typically used for Economic Dispatch where
# unit commitments are already fixed.

# ──► build_ptdf_model_with_quadratic_losses  [build_models.jl:895-916]
#
# Differences from the linear version:
#
#   STEP C2 is replaced by add_quadratic_current_loss_to_copperplate_balance!
#     [build_models.jl:340-359]
#     set_normalized_coefficient(CopperPlateBalanceConstraint[t], loss_var[t], 1)
#     (RHS is NOT shifted – it stays 0.  The quadratic constraint below
#      defines the loss variable's value instead of a fixed estimate.)
#     After: Σ injection[i,t] + loss_var[t] = 0
#
#   STEP C3 is replaced by add_current_loss_constraint_quadratic_approximation!
#     [build_models.jl:462-513]
#     Adds quadratic constraint (LineLossConstraintApproximation):
#     loss_var[t] = -Σ_k R[k] * (Σ_j V_line[k,t]/V_bus[j,t] * PTDF[k,j] * inj[j,t])²
#     where:
#       R[k]          = series resistance of branch k  (get_RX_vector, utils.jl:274-286)
#       V_line[k,t]   = max voltage at branch k endpoints  (utils.jl:242-254)
#       V_bus[j,t]    = bus voltage magnitude from prior AC PF  (utils.jl:201-209)
#       PTDF[k,j]     = power transfer distribution factor
#     Negative sign: losses consume power → reduce net injection available
#
# Note: STEP D (FND correction to branch constraints) is present in the code
# but commented out in this function because the quadratic constraint already
# captures the loss effect directly in the balance.

model_quad = build_ptdf_model_with_quadratic_losses(
    sys_uc,
    res_lossless;                          # previous result for voltage magnitudes
    device_models = DEFAULT_ED_MODELS,    # ED formulation (no commitment variables)
    optimizer = optimizer_with_attributes(Ipopt.Optimizer),
)
solve!(model_quad)
res_quad = OptimizationProblemResults(model_quad)

println("=== Quadratic loss model solved ===")
println("Objective: ", res_quad.optimizer_stats[1, "objective_value"])

# Inspect the new variable and constraint added to the JuMP model:
#   model_quad.internal.container.variables[
#       PSI.VariableKey{LineLossTotalApproximation, PSY.System}("")]
#   model_quad.internal.container.constraints[
#       InfrastructureSystems.Optimization.ConstraintKey{LineLossConstraintApproximation, PSY.System}("")]

# =============================================================================
# PART 4 – FLOW CANCELLING (transmission expansion, lossless)
# =============================================================================
# Flow cancelling enables an investment model where candidate lines can be
# "connected" (binary z_k = 1) or ignored (z_k = 0).  When z_k = 0 the
# candidate line still appears in the PTDF matrix (it was there at build!
# time), so its flow would normally affect other branches.  The flow-cancelling
# variables v_k[t] ≈ z_k * flow[k,t] are used to subtract out those spurious
# flows when the line is not built.
#
# Reference: Caramanis et al., IEEE PES 2016, eq. (3) and (18).

# Build the 5-bus investment system: existing generators + candidate generators
# + existing lines + candidate lines (parallel topology via intermediate buses)
sys_inv = build_matpower_5bus_with_updated_lines()
transform_single_time_series!(sys_inv, Hour(2), Hour(2))
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
model_fc = build_model_with_flow_canceling_terms(sys_inv)
solve!(model_fc)

res_fc = OptimizationProblemResults(model_fc)
println("=== Flow-cancelling model (lossless) solved ===")
println("Objective: ", JuMP.objective_value(model_fc.internal.container.JuMPmodel))

# Read investment decisions
inv_lines = read_variable(res_fc, PSI.VariableKey{BranchInvestmentVariable, Line}(""))
inv_gens  = read_variable(res_fc, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions:\n", inv_lines)
println("Generation investment decisions:\n", inv_gens)

# =============================================================================
# PART 5 – FLOW CANCELLING WITH QUADRATIC LOSS APPROXIMATION
# =============================================================================
# Combines the flow-cancelling investment model (Part 4) with the quadratic
# loss formulation (Part 3).  Requires an NLP-capable solver because of the
# quadratic loss constraints.
#
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
#   loss_var[t] = -Σ_k R[k] * (Σ_j PTDF[k,j] * injection[j,t])²
#   Flat-voltage version of the quadratic constraint (V ≡ 1 p.u.):
#   simpler than the voltage-scaled variant in build_models.jl Part 3 but
#   still physically meaningful for systems near nominal voltage.

model_fc_quad = build_model_with_flow_canceling_and_quadratic_losses(
    sys_inv;
    optimizer = optimizer_with_attributes(Ipopt.Optimizer),
)
solve!(model_fc_quad)

res_fc_quad = OptimizationProblemResults(model_fc_quad)
println("=== Flow-cancelling + quadratic losses model solved ===")
println("Objective: ", JuMP.objective_value(model_fc_quad.internal.container.JuMPmodel))

inv_lines_quad = read_variable(res_fc_quad, PSI.VariableKey{BranchInvestmentVariable, Line}(""))
inv_gens_quad  = read_variable(res_fc_quad, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
println("Line investment decisions (with losses):\n", inv_lines_quad)
println("Generation investment decisions (with losses):\n", inv_gens_quad)

# =============================================================================
# WHERE TO GO NEXT
# =============================================================================
# Use the function-name references above as entry points into the source.
# The most instructive functions to read in order are:
#
#  LOSSES
#  ├── build_ptdf_model_without_losses          build_models.jl:201
#  ├── add_current_loss_variables!              build_models.jl:234
#  ├── add_current_loss_to_copperplate_balance! build_models.jl:286
#  ├── add_current_loss_constraint_approximation!            :383
#  ├── add_current_loss_constraint_quadratic_approximation!  :462
#  ├── update_transmission_constraints_with_losses!          :722
#  └── run_iterative_linear_loss_model          run_models.jl:170
#
#  FLOW CANCELLING
#  ├── build_model_with_flow_canceling_terms    FlowCancelling/build_models.jl:497
#  ├── add_branch_investment_variables!                      :158
#  ├── add_branch_cancelling_flow_variables!                 :187
#  ├── add_bigM_linking_constraints!                         :222
#  ├── add_shift_terms_to_existing_line_constraints!         :266
#  ├── add_shift_terms_to_candidate_line_constraints!        :312
#  └── build_model_with_flow_canceling_and_quadratic_losses  :666
#
#  UTILITY / POST-PROCESSING
#  ├── get_bus_loss_factors                     utils.jl:73
#  ├── get_total_AC_loss                        utils.jl:418
#  ├── get_fictitious_nodal_demand_by_loss      utils.jl:806
#  └── get_RX_vector                            utils.jl:274
