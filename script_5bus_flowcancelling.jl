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
import PowerNetworkMatrices

include("Systems/5bus/build_5bus.jl")
include("SiennaScripts/FlowCancelling/build_models.jl")


sys = build_matpower_5bus_with_updated_lines()
transform_single_time_series!(sys, Hour(2), Hour(2))
set_available!(get_component(PhaseShiftingTransformer, sys, "bus-3-bus-4-i_5"), false)
candidate_gens = candidate_projects_data(sys)
for gen in candidate_gens
    add_component!(sys, gen)
end

add_candidate_line_data_without_parallel!(sys)

model = build_model_with_flow_canceling_terms(sys)
solve!(model)

res = OptimizationProblemResults(model)
obj_fun = JuMP.objective_function(model.internal.container.JuMPmodel)
inv_g = read_variable(res, PSI.VariableKey{GenerationInvestmentVariable, ThermalStandard}(""))
inv_l = read_variable(res, PSI.VariableKey{BranchInvestmentVariable, Line}(""))