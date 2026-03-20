using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using Ipopt
using Xpress
using Dates
using JuMP
using PowerFlows
using CSV
using JLD2
using Logging
using InfrastructureSystems
import PowerNetworkMatrices

const PSY = PowerSystems
const PSI = PowerSimulations

#this_path = @__DIR__
#CATS_path = joinpath(this_path, "CATS-CaliforniaTestSystem", "Sienna")
#include(joinpath(CATS_path, "build_CATS.jl"))
#### Run this once to convert CSV load data to JLD2 format ####
#include(joinpath(CATS_path, "convert_load_csv_to_jld2.jl"))
#### Run this to build and save the system ####
#system = build_CATS_system(first_order = true)
#to_json(system, "CATS_saved_reduced_sys.json")

# Load the saved system #
system = System("CATS_saved_reduced_sys.json")
transform_single_time_series!(
    system,
    Hour(2),  # horizon
    Hour(2),   # interval
)
ptdf = PTDF(system)
UC_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicUnitCommitment,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

ED_PTDF_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicDispatch,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

ED_MODELS = Dict(
    Line => StaticBranchUnbounded,
    Transformer2W => StaticBranchUnbounded,
    ThermalStandard => ThermalBasicDispatch,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    SynchronousCondenser => SynchronousCondenserBasicDispatch,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

include("SiennaScripts/build_models.jl")
include("SiennaScripts/utils.jl")
include("SiennaScripts/run_models.jl")
include("SiennaScripts/build_simulations.jl")
include("SiennaScripts/run_simulations.jl")
include("SiennaScripts/add_hvdc.jl")

add_internal_hvdc!(system)
solver_xpress = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.01, "MAXTIME" => 60*60)
solver_ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)

### Debug code ###

sim = build_uc_ed_simulation_with_acopf(
    system,
    system;
    uc_models = UC_MODELS,
    ed_models = ED_MODELS,
    ptdf_uc = ptdf,
    ptdf_ed = ptdf,
    uc_optimizer = solver_xpress,
    ed_optimizer = solver_ipopt,
)

execute!(sim)

uc_model = sim.models.decision_models[1]
uc = uc_model.internal.container
optimize!(uc.JuMPmodel)
injection_old_uc = deepcopy(
    JuMP.value.(
        uc.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
            ActivePowerBalance,
            ACBus,
        }(
            "",
        )]
    ).data,
)

sim_res = SimulationResults(sim)
res_old_uc = get_decision_problem_results(sim_res, "UC")
res_old_ed = get_decision_problem_results(sim_res, "ED")

v_mag_pf = read_realized_variable(res_old_ed, "PowerFlowVoltageMagnitude__ACBus"; table_format = TableFormat.WIDE)

loss_factors = get_bus_loss_factors(res_old_ed; slack_number = "1951")      # ∂Loss/∂P at each bus
total_loss_est = get_total_AC_loss_CATS(res_old_ed)       # Total AC losses in MW

sim_new = build_uc_ed_simulation_with_acopf(
    system,
    system;
    uc_models = UC_MODELS,
    ed_models = ED_MODELS,
    ptdf_uc = ptdf,
    ptdf_ed = ptdf,
    uc_optimizer = solver_xpress,
    ed_optimizer = solver_ipopt,
)

uc_model_new = sim_new.models.decision_models[1]

update_copperplate_loss_approximation!(
        uc_model,
        loss_factors,
        total_loss_est,
        injection_old_uc,
    )

execute!(sim_new)
