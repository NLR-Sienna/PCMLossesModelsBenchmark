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
import PowerNetworkMatrices

this_path = @__DIR__
CATS_path = joinpath(this_path, "CATS-CaliforniaTestSystem", "Sienna")
include(joinpath(CATS_path, "build_CATS.jl"))
#### Run this once to convert CSV load data to JLD2 format ####
# include(joinpath(CATS_path, "convert_load_csv_to_jld2.jl"))
#### Run this to build and save the system ####
# system = build_CATS_system(first_order = true)
# to_json(system, "CATS_saved_sys.json")

# Load the saved system #
system = System("CATS_saved_sys.json")

transform_single_time_series!(
    system,
    Hour(1),  # horizon
    Hour(1),   # interval
)
ptdf = VirtualPTDF(system;
    tol=0.0001,
    max_cache_size=10000,
    # radial_network_reduction = RadialNetworkReduction(PNM.IncidenceMatrix(sys)), #Jose's idea
)

template_uc = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=true, duals=[CopperPlateBalanceConstraint]))
 
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, Line, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)
 
solver_xpress = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.01, "MAXTIME" => 60*60)
 
problem_uc = DecisionModel(
    template_uc,
    system;
    optimizer=solver_xpress,
    optimizer_solve_log_print=true,
    calculate_conflict=true,
    name="UC"
)

build!(problem_uc, output_dir=mktempdir())

solve!(problem_uc)

results_uc_1 = OptimizationProblemResults(problem_uc)

slack_up = read_variable(results_uc_1, "SystemBalanceSlackUp__System")
slack_dn = read_variable(results_uc_1, "SystemBalanceSlackDown__System")

##### Do ACOPF Sequence #####

# UC #
template_uc = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=true, duals=[CopperPlateBalanceConstraint]))
 
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, Line, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)
solver_xpress = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.01, "MAXTIME" => 60*60)
 
problem_uc = DecisionModel(
    template_uc,
    system;
    optimizer=solver_xpress,
    optimizer_solve_log_print=true,
    calculate_conflict=true,
    name="UC"
)

# ED ACOPF #
network_model_ed = NetworkModel(ACPPowerModel; use_slacks=true, power_flow_evaluation = PowerFlows.ACPowerFlow(; calculate_loss_factors=true, calculate_voltage_stability_factors = true))
template_ed = ProblemTemplate(network_model_ed)
set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_device_model!(template_ed, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_ed, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_ed, PowerLoad, StaticPowerLoad)
set_device_model!(template_ed, Line, StaticBranchUnbounded)
set_device_model!(template_ed, Transformer2W, StaticBranchUnbounded)
set_device_model!(template_ed, SynchronousCondenser, SynchronousCondenserBasicDispatch)

solver_ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)
problem_ed = DecisionModel(
    template_ed, 
    system; 
    optimizer = solver_ipopt,
    optimizer_solve_log_print = true, 
    name = "ED"
)

models = SimulationModels(;
    decision_models = [
        problem_uc,
        problem_ed,
    ],
)

sequence = SimulationSequence(;
    models = models,
    feedforwards = Dict(
        "ED" => [
            SemiContinuousFeedforward(;
                component_type = ThermalStandard,
                source = OnVariable,
                affected_values = [ActivePowerVariable],
            ),
        ],
    ),
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(;
    name = "no_cache",
    steps = 2,
    models = models,
    sequence = sequence,
    simulation_folder = mktempdir(),
)

build!(sim)
execute!(sim)

sim_res = SimulationResults(sim)
results_acopf = get_decision_problem_results(sim_res, "ED")


v_mag_pf = read_realized_variable(results_acopf, "PowerFlowVoltageMagnitude__ACBus"; table_format = TableFormat.WIDE)
loss_factors = read_realized_variable(results_acopf, "PowerFlowLossFactors__ACBus"; table_format = TableFormat.WIDE)
stability_factors = read_realized_variable(results_acopf, "PowerFlowVoltageStabilityFactors__ACBus"; table_format = TableFormat.WIDE)
loss_factors[!, 1000:1020]
p_slack_up = read_realized_variable(results_acopf, "SystemBalanceSlackDown__ACBus__P"; table_format = TableFormat.WIDE)
p_slack_up_tot = sum(eachcol(p_slack_up[!, 2:end]))

p_slack_dn = read_realized_variable(results_acopf, "SystemBalanceSlackDown__ACBus__P"; table_format = TableFormat.WIDE)
p_slack_dn_tot = sum(eachcol(p_slack_dn[!, 2:end]))

q_slack_up = read_realized_variable(results_acopf, "SystemBalanceSlackUp__ACBus__Q"; table_format = TableFormat.WIDE)
q_slack_up_tot = sum(eachcol(q_slack_up[!, 2:end]))
q_slack_dn = read_realized_variable(results_acopf, "SystemBalanceSlackDown__ACBus__Q"; table_format = TableFormat.WIDE)
q_slack_dn_tot = sum(eachcol(q_slack_dn[!, 2:end]))

q_sl = q_slack_up[!, 2:end] .- q_slack_dn[!, 2:end]
q_sl_tot = sum(eachcol(q_sl))       