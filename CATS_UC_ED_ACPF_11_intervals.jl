
using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using Ipopt
using Gurobi
using Dates
using JuMP
using PowerFlows
import PowerNetworkMatrices: VirtualPTDF


#=
using Pkg
Pkg.status("PowerSimulations") 
Pkg.status("PowerFlows.jl") 
Pkg.develop(url="https://github.com/NREL-Sienna/PowerSimulations.jl", rev="rb/loss_factors_variable")
Pkg.add(url="https://github.com/NREL-Sienna/PowerFlows.jl.git", rev="lin-solve-cache")
=#

const PSI = PowerSimulations

#import Xpress

import HSL_jll # works only after getting the HSL license and the files for HSL_jll
using HSL # works only after getting the HSL license and the files for HSL_jll (]dev HSL_jll first)

mip_gap = 0.5

include("CATS-CaliforniaTestSystem/build-system/parse-matpower.jl");
#system = System("CATS-CaliforniaTestSystem/MATPOWER/CaliforniaTestSystem.m");
system = System("CATS-CaliforniaTestSystem/MATPOWER/system_condensers_removed_cutoff75.m");
include("CATS-CaliforniaTestSystem/build-system/replace_gens.jl");
include("CATS-CaliforniaTestSystem/build-system/define_time_series.jl");

system_ed = deepcopy(system)

# for multiple time steps:
transform_single_time_series!(
    system,
    Hour(11),  # horizon
    Hour(11),   # interval
);

transform_single_time_series!(
    system_ed,
    Hour(1),  # horizon
    Hour(1),  # interval
);

ptdf = VirtualPTDF(system;
    tol=0.0001,
    max_cache_size=10000,
    # radial_network_reduction = RadialNetworkReduction(PNM.IncidenceMatrix(sys)), #Jose's idea
)

network_model_uc = NetworkModel(PTDFPowerModel; PTDF_matrix=ptdf)
network_model_ed = NetworkModel(ACPPowerModel; use_slacks=true, power_flow_evaluation=PowerFlows.ACPowerFlow(; calculate_loss_factors=true))
#network_model_ed = NetworkModel(ACPPowerModel; use_slacks=true)

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

solver_gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer)

solver_ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    "print_level" => 5,
    "hsllib" => HSL_jll.libhsl_path, # uncomment after getting the HSL license and the files for HSL_jll
    "linear_solver" => "ma57", # uncomment after getting the HSL license and the files for HSL_jll
    "tol" => 1e-3,
    "acceptable_tol" => 1e-3,
)


problem_uc = DecisionModel(
    template_uc,
    system;
    optimizer=solver_gurobi,
    optimizer_solve_log_print=true,
    name="UC"
)

problem_ed = DecisionModel(
    template_ed,
    system_ed;
    optimizer=solver_ipopt,
    optimizer_solve_log_print=true,
    name="ED"
)

models = SimulationModels(;
    decision_models=[
        problem_uc,
        problem_ed,
    ],
)

sequence = SimulationSequence(;
    models=models,
    feedforwards=Dict(
        "ED" => [
            SemiContinuousFeedforward(;
                component_type=ThermalStandard,
                source=OnVariable,
                affected_values=[ActivePowerVariable],
            ),
        ],
    ),
    ini_cond_chronology=InterProblemChronology(),
)

sim = Simulation(;
    name="no_cache",
    steps=2,
    models=models,
    sequence=sequence,
    simulation_folder=mktempdir(),
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

# Open a file to capture output
open("ipopt_log.txt", "w") do io
    redirect_stdout(io)  # Redirect standard output
    redirect_stderr(io)  # Redirect error output (failure messages)

    try
        execute_out = execute!(sim; exports=exports, in_memory=true)
        @assert execute_out == PSI.RunStatus.SUCCESSFULLY_FINALIZED
    catch err
        println("Solver failed: ", err)
    end
end


results = SimulationResults(sim);
uc_results = get_decision_problem_results(results, "UC")
ed_results = get_decision_problem_results(results, "ED")

vduc = read_variables(uc_results)
on = vduc["OnVariable__ThermalStandard"]

vd = read_variables(ed_results)

ad = read_aux_variables(ed_results)

# reading the loss factors from ED results:
lf_res = ad["PowerFlowLossFactors__ACBus"]

slackdn_q = read_variable(ed_results, "SystemBalanceSlackDown__ACBus__Q")
slackup_q = read_variable(ed_results, "SystemBalanceSlackUp__ACBus__Q")
slackdn_p = read_variable(ed_results, "SystemBalanceSlackDown__ACBus__P")
slackup_p = read_variable(ed_results, "SystemBalanceSlackUp__ACBus__P")

# Check buses with slacks
slack_down_q = Set()
for (k, v) in slackdn_q
    for n in names(slackdn_q[k])
        try
            #ss=round(sum(slackdn_q[k][!,n]),digits=2)
            ss = round(sum(slackdn_q[k][:, n]), digits=2)
            if abs(ss) > 0.00001
                println("k,", k, ",n,", n, ",", ss)
                push!(slack_down_q, n)
            end
        catch e
        end
    end
end
println(slack_down_q)

slack_up_q = Set()
for (k, v) in slackup_q
    for n in names(slackup_q[k])
        try
            #ss=round(sum(slackup_q[k][!,n]),digits=2)
            ss = round(sum(slackup_q[k][:, n]), digits=2)
            if abs(ss) > 0.00001
                println("k,", k, ",n,", n, ",", ss)
                push!(slack_up_q, n)
            end
        catch e
        end
    end
end
println(slack_up_q)

slack_down_p = Set()
for (k, v) in slackdn_p
    for n in names(slackdn_p[k])
        try
            #ss=round(sum(slackdn_p[k][!,n]),digits=2)
            ss = round(sum(slackdn_p[k][:, n]), digits=2)
            if abs(ss) > 0.00001
                println("k,", k, ",n,", n, ",", ss)
                push!(slack_down_p, n)
            end
        catch e
        end
    end
end
println(slack_down_p)

slack_up_p = Set()
for (k, v) in slackup_p
    for n in names(slackup_p[k])
        try
            #ss=round(sum(slackup_p[k][!,n]),digits=2)
            ss = round(sum(slackup_p[k][:, n]), digits=2)
            if abs(ss) > 0.00001
                println("k,", k, ",n,", n, ",", ss)
                push!(slack_up_p, n)
            end
        catch e
        end
    end
end
println(slack_up_p)