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
using CSV
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
PowerSystemCaseBuilder.clear_all_serialized_systems()
system_name = "matpower_case5_sys"
sys = build_system(MatpowerTestSystems, system_name)
transform_single_time_series!(sys, Hour(1), Hour(1))
sys_ts = System("modified_RTS_GMLC_DA_sys_noForecast.json")
load_sys = ["bus2", "bus3", "bus4"]
load_sys_ts = ["Baker", "Bacon", "Bajer"]
for (l, lt) in zip(load_sys, load_sys_ts)
    load = get_component(StaticLoad, sys, l)
    load_ts = get_component(StaticLoad, sys_ts, lt)
    ts_array = get_time_series_array(SingleTimeSeries, load_ts, "max_active_power"; ignore_scaling_factors = true)
    values(ts_array)[1] = 1.0 # Set load to 1.0 p.u. for first time step
    new_ts = SingleTimeSeries(; name = "max_active_power", data = ts_array, scaling_factor_multiplier=get_max_active_power)
    add_time_series!(sys, load, new_ts)
end
update_generation_costs!(sys)
update_line_ratings!(sys, 1.0)

output_dir = "./5bus_PCM"

optimizer =  optimizer_with_attributes(Gurobi.Optimizer)
const DEFAULT_UC_MODELS = Dict(
    Line => StaticBranchBounds,
    PhaseShiftingTransformer => StaticBranch,
    ThermalStandard => ThermalDispatchNoMin,
    PowerLoad => StaticPowerLoad,
)
device_models = DEFAULT_UC_MODELS
template = ProblemTemplate(
            NetworkModel(
                ACPPowerModel;
                use_slacks = true,
            ),)
for (device_type, formulation) in device_models
    set_device_model!(template, device_type, formulation)
end

model = DecisionModel(
    template,
    sys;
    optimizer = optimizer,
    name = "ED",
    store_variable_names = true,
)
models = SimulationModels(
           decision_models=[model]);
sequence =
           SimulationSequence(models=models, ini_cond_chronology=InterProblemChronology());
sim = Simulation(
           name="5bus_PCM_AC",
           steps = 8760,
           models = models,
           sequence = sequence,
           simulation_folder=output_dir);
build!(sim)
execute!(sim)

function save_to_csv(cont, model_name, save_dir=".")
           for (name, df) in cont
               CSV.write(joinpath(save_dir, "$(name)_$(model_name).csv"), df)
           end
       end
function export_results_csv(results, variable, stage, path)
           variables = variable
           aux_variables = PSI.read_realized_aux_variables(results)
           parameters = PSI.read_realized_parameters(results)
           duals = PSI.read_realized_duals(results)
           expressions = PSI.read_realized_expressions(results)

           save_to_csv(variables, stage, path)
           save_to_csv(parameters, stage, path)
           save_to_csv(duals, stage, path)
           save_to_csv(aux_variables, stage, path)
           save_to_csv(expressions, stage, path)
           return
       end

model = get_simulation_model(sim, :ED);
results = SimulationResults(sim; ignore_status=true);
results_ed = get_decision_problem_results(results, "ED")
set_system!(results_ed, sys)
variables = PSI.read_realized_variables(results_ed)
duals = PSI.read_realized_duals(results_ed)
expressions = PSI.read_realized_expressions(results_ed)
parameters = PSI.read_realized_parameters(results_ed)
export_results_csv(results_ed, variables, "ED", joinpath(results.path, "results"))