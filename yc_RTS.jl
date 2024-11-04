using Pkg
Pkg.activate(".")
using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using Ipopt
using PowerFlows
using TimeSeries
using HiGHS
using Logging
using Revise
using Xpress
using JuMP
using Dates
using CSV
using DataFrames
using PowerGraphics
#in pkg dev PowerSimulations

const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
const PSI = PowerSimulations

#sys = build_system(PSISystems, "c_sys5_pjm")
PSB.clear_all_serialized_systems() 

c_sys5_pjm_da = PSB.build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization = true)
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(24), Hour(24))
c_sys5_pjm_rt = PSB.build_system(PSISystems, "modified_RTS_GMLC_RT_sys"; skip_serialization = true)
#c_sys5_pjm_rt = PSB.build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization = true)
#for x in get_components(ACBus, c_sys5_pjm_rt)
#    set_voltage_limits!(x, (min = 0.9, max = 1.1))
#end
PSY.transform_single_time_series!(c_sys5_pjm_rt, Hour(1), Hour(1))

#sys_DA = build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization = true)
#sys_RT = build_system(PSISystems, "modified_RTS_GMLC_RT_sys"; skip_serialization = true)
#=
template_uc = template_unit_commitment()
pop!(template_uc.branches, :TwoTerminalHVDCLine)
set_network_model!(template_uc,    NetworkModel(PTDFPowerModel;),)
#set_network_model!(template_uc,    NetworkModel(DCPLLPowerModel;),)
set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
=#

template_uc = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=true))
set_device_model!(template_uc, Line, StaticBranchUnbounded)
set_device_model!(template_uc, Transformer2W, StaticBranchUnbounded)
set_device_model!(template_uc, TapTransformer, StaticBranchUnbounded)
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
#set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, RenewableNonDispatch, FixedOutput)

template_ed = deepcopy(template_uc)
set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_network_model!(template_ed, NetworkModel( ACPPowerModel,use_slacks = true ),)

#milp_optimizer = optimizer_with_attributes(Xpress.Optimizer)
milp_optimizer=optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.02)
nlp_optimizer = optimizer_with_attributes(Ipopt.Optimizer)

#= setting must run can address slackdn_q
mustrun_gen_set=["307_CT_2","207_CT_2","207_CT_1","307_CT_1"]
for m in mustrun_gen_set
    g=get_component(ThermalStandard, c_sys5_pjm_da, m)
    set_must_run!(g, true)
    g=get_component(ThermalStandard, c_sys5_pjm_rt, m)
    set_must_run!(g, true)
end
=#

models = SimulationModels(
    decision_models=[
        DecisionModel(
            template_uc,
            c_sys5_pjm_da,
            optimizer=milp_optimizer,
            name="UC",
#            system_to_file=true,
            store_variable_names=true,
        ),
        DecisionModel(
            template_ed,
            c_sys5_pjm_rt,
            optimizer=nlp_optimizer,
            name="ED",
#            system_to_file=true,
            initialize_model=false,
            store_variable_names=true,
        ),
    ],
)

sequence = SimulationSequence(
    models=models,
    feedforwards=Dict(
        "ED" => [SemiContinuousFeedforward(component_type=ThermalStandard,  source=OnVariable,
                affected_values=[ActivePowerVariable, ReactivePowerVariable],
                #affected_values=[ActivePowerVariable],
            ),
        ],
    ),
    ini_cond_chronology=InterProblemChronology(),
)

sim = Simulation(
    name="RTS_sim",
    steps=1,
    models=models,
    sequence=sequence,
    simulation_folder=".",
)

build!(sim; console_level=Logging.Info)
execute!(sim)

res_sim_yc = SimulationResults(sim)
results_uc_yc = get_decision_problem_results(res_sim_yc, "UC")
results_ed_yc = get_decision_problem_results(res_sim_yc, "ED")

plot_fuel(results_uc_yc; stair = true)
plot_fuel(results_ed_yc; stair=true)

vmag_yc= read_variable(results_ed_yc, "VoltageMagnitude__ACBus")

th_on_yc = read_variable(results_uc_yc, "OnVariable__ThermalStandard")
th_power_yc = read_variable(results_ed_yc, "ActivePowerVariable__ThermalStandard")
th_rpower_yc = read_variable(results_ed_yc, "ReactivePowerVariable__ThermalStandard")
load_param_yc = read_parameter(results_ed_yc, "ActivePowerTimeSeriesParameter__PowerLoad")
load_param_re_yc = read_parameter(results_ed_yc, "ReactivePowerTimeSeriesParameter__PowerLoad")

slackdn_q=read_variable(results_ed_yc, "SystemBalanceSlackDown__ACBus__Q")
slackup_q=read_variable(results_ed_yc, "SystemBalanceSlackUp__ACBus__Q")
slackdn_p=read_variable(results_ed_yc, "SystemBalanceSlackDown__ACBus__P")
slackup_p=read_variable(results_ed_yc, "SystemBalanceSlackUp__ACBus__P")

# Check buses with slacks
s=Set()
for (k,v) in slackdn_q
    for n in names(slackdn_q[k]) 
        try 
            ss=round(sum(slackdn_q[k][!,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(s,n)
            end 
        catch e end  
    end 
end
println(s)

for (k,v) in slackup_q
    for n in names(slackup_q[k]) 
        try 
            ss=round(sum(slackup_q[k][!,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(s,n)
            end 
        catch e end  
    end 
end
println(s)

for (k,v) in slackdn_p
    for n in names(slackdn_p[k]) 
        try 
            ss=round(sum(slackdn_p[k][!,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
            end 
        catch e end  
    end 
end

for (k,v) in slackup_p
    for n in names(slackup_p[k]) 
        try 
            ss=round(sum(slackup_p[k][!,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
            end 
        catch e end  
    end 
end

for b in PSY.get_components(PSY.Generator, c_sys5_pjm_da)
    if string(get_number(get_bus(b))) in s
        println(get_number(get_bus(b)),",",get_name(b))
    end        
end 

#=
308,308_RTPV_1
307,307_CT_2
207,207_CT_2
207,207_CT_1
307,307_CT_1
=#

