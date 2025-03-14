using AppleAccelerate
using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
# using HiGHS # solver
using Ipopt
using Dates
using JuMP
using PowerFlows
import PowerNetworkMatrices: VirtualPTDF

const PSI = PowerSimulations

import Xpress

import HSL_jll # works only after getting the HSL license and the files for HSL_jll
using HSL # works only after getting the HSL license and the files for HSL_jll (]dev HSL_jll first)

mip_gap = 0.5

PROJECT_ROOT = realpath(".")  # TODO set the path to the project directory of CATS-CaliforniaTestSystem
build_system_path = joinpath(PROJECT_ROOT, "build-system", "build_from_matpower.jl")
include(build_system_path)
display(system)

system_ed = deepcopy(system)

# # for multiple time step horizon (slow):
# transform_single_time_series!(
#            system,
#            Dates.Hour(48),  # horizon: 48 hr ahead 
#            Dates.Hour(24),   # interval 
#        );

# for one time step only:
transform_single_time_series!(
           system,
           Dates.Hour(1),  # horizon: 48 hr ahead 
           Dates.Hour(1),   # interval 
       );

transform_single_time_series!(
        system_ed,
        Dates.Hour(1),  # horizon: 1 hr ahead
        Dates.Hour(1),  # interval 
    );

ptdf = VirtualPTDF(system;
        tol = .0001,
        max_cache_size = 10000,
        # radial_network_reduction = RadialNetworkReduction(PNM.IncidenceMatrix(sys)), #Jose's idea
        )

network_model_uc = NetworkModel(PTDFPowerModel; PTDF_matrix=ptdf)
network_model_ed = NetworkModel(ACPPowerModel; use_slacks=true, power_flow_evaluation=PowerFlows.ACPowerFlow(;calculate_loss_factors=true))

template_uc = ProblemTemplate(network_model_uc)
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, Line, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)


template_ed = ProblemTemplate(network_model_ed)
set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_device_model!(template_ed, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_ed, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_ed, PowerLoad, StaticPowerLoad)
set_device_model!(template_ed, Line, StaticBranchUnbounded)
set_device_model!(template_ed, Transformer2W, StaticBranchUnbounded)

# solver_highs = optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => mip_gap) #, "presolve" => "off" ) 

solver_xpress = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.01)

solver_ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    "print_level" => 3,
    "hsllib" => HSL_jll.libhsl_path, # uncomment after getting the HSL license and the files for HSL_jll
    "linear_solver" => "ma57", # uncomment after getting the HSL license and the files for HSL_jll
    "tol" => 1e-6,
    "acceptable_tol" => 1e-3,
)


problem_uc = DecisionModel(
    template_uc, 
    system; 
    optimizer = solver_xpress, 
    optimizer_solve_log_print = true, 
    name = "UC"
)

problem_ed = DecisionModel(
    template_ed, 
    system_ed; 
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

build_out = build!(sim)
@assert build_out == PSI.SimulationBuildStatus.BUILT

exports = Dict(
    "models" => [
        Dict(
            "name" => "UC",
            "store_all_variables" => true,
            "store_all_parameters" => true,
            "store_all_duals" => true,
            "store_all_aux_variables" => true,
        ),
        Dict(
            "name" => "ED",
            "store_all_variables" => true,
            "store_all_parameters" => true,
            "store_all_duals" => true,
            "store_all_aux_variables" => true,
        ),
    ],
    "path" => mktempdir(),
    "optimizer_stats" => true,
)
execute_out = execute!(sim; exports = exports, in_memory = true)
@assert execute_out == PSI.RunStatus.SUCCESSFULLY_FINALIZED

results = SimulationResults(sim);
uc_results = get_decision_problem_results(results, "UC")
ed_results = get_decision_problem_results(results, "ED")

vd = read_variables(ed_results)

ad = read_aux_variables(ed_results)

# reading the loss factors from ED results:
lf_res=ad["PowerFlowLossFactors__ACBus"]