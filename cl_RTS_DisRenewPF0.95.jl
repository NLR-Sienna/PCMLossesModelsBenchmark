using Pkg
#Pkg.actiaddvate(".")
#Pkg.activate(".")
Pkg.activate(@__DIR__)
using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using Ipopt
using PowerFlows
using TimeSeries
using HiGHS
using Logging
using Revise
#using Xpress
using JuMP
using Dates
using CSV
using DataFrames
using PowerGraphics
#in pkg dev PowerSimulations

const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
const PSI = PowerSimulations

#clear_all_serialized_systems() 

c_sys5_pjm_da = PSB.build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization = true)
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(24), Hour(24))

c_sys5_pjm_rt = PSB.build_system(PSISystems, "modified_RTS_GMLC_RT_sys"; skip_serialization = true)

#=
for x in get_components(ACBus, c_sys5_pjm_rt)
    set_voltage_limits!(x, (min = 0.94, max = 1.05))
end
=#


# Add the upper bound of reactive power for renewable energy resources to solve the votlage problem
#=
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "309_WIND_1"), (min=0, max=4.6651))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "303_WIND_1"), (min=0, max=4.6651))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "317_WIND_1"), (min=0, max=4.6651))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "122_WIND_1"), (min=0, max=4.6651))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "103_PV_1"), (min=0, max=4.6651))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "104_PV_1"), (min=0, max=4.6651))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "324_PV_1"), (min=0, max=4.6651))
=#

Renewable = collect(get_components(RenewableDispatch, c_sys5_pjm_rt))
Num = length(Renewable)
for row in 1:Num
    set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, get_name(Renewable[row])), 0.95)
    R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,get_name(Renewable[row])))
    PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,get_name(Renewable[row])))
    Q_max = R * sin(acos(PF))
    Q_min = -R * sin(acos(PF))
    #Q_min = 0
    set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, get_name(Renewable[row])), (min=Q_min, max=Q_max))
end

#=
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "309_WIND_1"), 0.95)
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"309_WIND_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"309_WIND_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "309_WIND_1"), (min=Q_min, max=Q_max))
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "303_WIND_1"), 0.95)
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"303_WIND_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"303_WIND_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "303_WIND_1"), (min=Q_min, max=Q_max))
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "317_WIND_1"), 0.95)
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"317_WIND_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"317_WIND_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "317_WIND_1"), (min=Q_min, max=Q_max))
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "122_WIND_1"), 0.95)
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"122_WIND_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"122_WIND_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "122_WIND_1"), (min=Q_min, max=Q_max))
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "103_PV_1"), 0.95)
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"103_PV_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"103_PV_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "103_PV_1"), (min=Q_min, max=Q_max))
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "104_PV_1"), 0.95)
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"104_PV_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"104_PV_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "104_PV_1"), (min=Q_min, max=Q_max))
set_power_factor!(get_component(RenewableDispatch, c_sys5_pjm_rt, "324_PV_1"), 0.95)
#get_max_active_power(get_component(RenewableDispatch,c_sys5_pjm_rt,"324_PV_1"))
R = get_rating(get_component(RenewableDispatch,c_sys5_pjm_rt,"324_PV_1"))
PF = get_power_factor(get_component(RenewableDispatch,c_sys5_pjm_rt,"324_PV_1"))
Q_max = R * sin(acos(PF))
Q_min = -R * sin(acos(PF))
set_reactive_power_limits!(get_component(RenewableDispatch, c_sys5_pjm_rt, "324_PV_1"), (min=Q_min, max=Q_max))
#get_max_reactive_power(get_component(RenewableDispatch,c_sys5_pjm_rt,"324_PV_1"))
=#
PSY.transform_single_time_series!(c_sys5_pjm_rt, Hour(1), Hour(1))


template_uc = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=false))
set_device_model!(template_uc, Line, StaticBranchUnbounded)
set_device_model!(template_uc, Transformer2W, StaticBranchUnbounded)
set_device_model!(template_uc, TapTransformer, StaticBranchUnbounded)
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)



template_ed = deepcopy(template_uc)
set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_network_model!(template_ed, NetworkModel( ACPPowerModel,use_slacks = true ),)

#milp_optimizer=optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.02)
milp_optimizer=optimizer_with_attributes(HiGHS.Optimizer)
nlp_optimizer = optimizer_with_attributes(Ipopt.Optimizer)

#=
# setting must run can address slackup_q
mustrun_gen_set=["307_CT_2","207_CT_2","207_CT_1","307_CT_1"]
for m in mustrun_gen_set
    g=get_component(ThermalStandard, c_sys5_pjm_da, m)
    set_must_run!(g, true)
    g=get_component(ThermalStandard, c_sys5_pjm_rt, m)
    set_must_run!(g, true)
end
=#

#=
for b in PSY.get_components(PSY.RenewableDispatch, c_sys5_pjm_da)
    set_available!(b,false)    
end 
for b in PSY.get_components(PSY.RenewableDispatch, c_sys5_pjm_rt)
    set_available!(b,false)    
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

#=sim = Simulation(
    name="RTS_sim",
    steps=1,
    models=models,
    sequence=sequence,
    simulation_folder=".",
) =#

sim = Simulation(
    name="RTS_sim",
    steps=1,
    models=models,
    sequence=sequence,
    simulation_folder=mktempdir(),
)


build!(sim; console_level=Logging.Info)
#build!(sim)
execute!(sim)


## Output results
# Retriving results
res_sim_yc = SimulationResults(sim)
results_uc_yc = get_decision_problem_results(res_sim_yc, "UC")
results_ed_yc = get_decision_problem_results(res_sim_yc, "ED")

#plot_fuel(results_uc_yc; stair = true)
#plot_fuel(results_ed_yc; stair=true)

plot_fuel(results_uc_yc; stair = true, format = png)
plot_fuel(results_ed_yc; stair=true, format = png)

vmag_yc = read_variable(results_ed_yc, "VoltageMagnitude__ACBus")
flag_up = 1 #if the nodal volatge is greater than 1.05, flag_up=0, otherwise, flag_up=1
flag_down = 1 #if the nodal volatge is lower than 0.95, flag_down=0, otherwise, flag_down=1
V_up = Set()
V_down = Set()
for (datetime, df) in vmag_yc
    for row in 1:nrow(df)
        for col in 2:ncol(df)
            vv1 = round(float(df[row, col]),digits=2)
            vv2 = float(df[row, col])
            #vv = float(df[row, col])
            #println(vv)
            if vv1 > 1.05
                flag_up = 0
                vv_coll = [vv1, col]
                #println(vv1)
                push!(V_up,vv_coll')
            end
            if vv2 < 0.949
                #println(vv2)
                vv_col = [vv2, col]
                flag_down = 0
                push!(V_down,vv_col')
            end
        end
    end
end
println(flag_up)
println(flag_down)

#= Output the V_down to Excel
a = collecr(V_down)
matrix = reshape(a,392,1)
XLSX.openxlsx("output.xlsx", mode="w") do xf
    # Write matrix data to the default sheet "Sheet1"
    sheet = xf[1]  # Access the first sheet by index
    for i in 1:size(matrix, 1)
        new_matrix = Matrix(matrix[i])
        for j in 1:size(new_matrix, 2)
            sheet[i, j] = new_matrix[1, j]
        end
    end
end
=#


th_on_yc = read_variable(results_uc_yc, "OnVariable__ThermalStandard")
th_power_yc = read_variable(results_ed_yc, "ActivePowerVariable__ThermalStandard")
th_rpower_yc = read_variable(results_ed_yc, "ReactivePowerVariable__ThermalStandard")
re_power_yc = read_variable(results_ed_yc, "ActivePowerVariable__RenewableDispatch")
re_rpower_yc = read_variable(results_ed_yc, "ReactivePowerVariable__RenewableDispatch")
load_param_yc = read_parameter(results_ed_yc, "ActivePowerTimeSeriesParameter__PowerLoad")
load_param_re_yc = read_parameter(results_ed_yc, "ReactivePowerTimeSeriesParameter__PowerLoad")

slackdn_q=read_variable(results_ed_yc, "SystemBalanceSlackDown__ACBus__Q")
slackup_q=read_variable(results_ed_yc, "SystemBalanceSlackUp__ACBus__Q")
slackdn_p=read_variable(results_ed_yc, "SystemBalanceSlackDown__ACBus__P")
slackup_p=read_variable(results_ed_yc, "SystemBalanceSlackUp__ACBus__P")


# check the dispatched renewable energy resources during ED
Ra_output = Set()
#=
for (datetime, df) in re_power_yc
    for row in 1:nrow(df)
        for col in 2:ncol(df)
            ou = round(float(df[row, col]),digits=2)
            if ou > 0.1
                ou_coll = [ou, col]
                #println(vv1)
                push!(Ra_output,ou_coll')
            end
        end
    end
end
=#
for (k,o) in re_power_yc
    for n in names(re_power_yc[k]) 
        try 
            ss=round(sum(re_power_yc[k][:,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(Ra_output,n)
            end 
        catch e end  
    end 
end


Rr_output = Set()
#=
for (datetime, df) in re_rpower_yc
    for row in 1:nrow(df)
        for col in 2:ncol(df)
            rou = abs(round(float(df[row, col]),digits=2))
            if rou > 0.0001
                rou_coll = [rou, col]
                #println(vv1)
                push!(Rr_output,rou_coll')
            end
        end
    end
end
=#
for (k,o) in re_rpower_yc
    for n in names(re_rpower_yc[k]) 
        try 
            ss=round(sum(re_rpower_yc[k][:,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(Rr_output,n)
            end 
        catch e end  
    end 
end


# Check buses with slacks
slack_down_q=Set()
for (k,v) in slackdn_q
    for n in names(slackdn_q[k]) 
        try 
            #ss=round(sum(slackdn_q[k][!,n]),digits=2)
            ss=round(sum(slackdn_q[k][:,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(slack_down_q,n)
            end 
        catch e end  
    end 
end
println(slack_down_q)

slack_up_q=Set()
for (k,v) in slackup_q
    for n in names(slackup_q[k]) 
        try 
            #ss=round(sum(slackup_q[k][!,n]),digits=2)
            ss=round(sum(slackup_q[k][:,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(slack_up_q,n)
            end 
        catch e end  
    end 
end
println(slack_up_q)

slack_down_p=Set()
for (k,v) in slackdn_p
    for n in names(slackdn_p[k]) 
        try 
            #ss=round(sum(slackdn_p[k][!,n]),digits=2)
            ss=round(sum(slackdn_p[k][:,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(slack_down_p,n)
            end 
        catch e end  
    end 
end
println(slack_down_p)

slack_up_p=Set()
for (k,v) in slackup_p
    for n in names(slackup_p[k]) 
        try 
            #ss=round(sum(slackup_p[k][!,n]),digits=2)
            ss=round(sum(slackup_p[k][:,n]),digits=2)
            if abs(ss)>0.00001
             println("k,",k,",n,",n,",",ss) 
             push!(slack_up_p,n)
            end 
        catch e end  
    end 
end
println(slack_up_p)

# check the generators at the buses with slacks
flag=0
for b in PSY.get_components(PSY.Generator, c_sys5_pjm_da)
    if string(get_number(get_bus(b))) in slack_down_q
        println(get_number(get_bus(b)),",",get_name(b))
        flag=1
        println(flag)
    elseif string(get_number(get_bus(b))) in slack_up_q
        println(get_number(get_bus(b)),",",get_name(b)) 
        flag=2
        println(flag) 
    elseif string(get_number(get_bus(b))) in slack_down_p
        println(get_number(get_bus(b)),",",get_name(b)) 
        flag=3
        println(flag)
    elseif string(get_number(get_bus(b))) in slack_up_p
        println(get_number(get_bus(b)),",",get_name(b))  
        flag=4
        println(flag)              
    end        
end 

#=
308,308_RTPV_1
307,307_CT_2
207,207_CT_2
207,207_CT_1
307,307_CT_1
=#




#compare the load and output of the generators
indexx = DateTime("2020-01-01T00:00:00")
load_in = load_param_yc[indexx]
Thermal_output = th_power_yc[indexx]
Renewable_output = re_power_yc[indexx]
slack_dp = slackdn_p[indexx]
slack_up = slackup_p[indexx]
load_q = load_param_re_yc[indexx]
Thermal_routput = th_rpower_yc[indexx]
Renewable_routput = re_rpower_yc[indexx]
slack_dq = slackdn_q[indexx]
slack_uq = slackup_q[indexx]
h = 1
sum_load = sum(load_in[h,2:end])
sum_thermal = sum(Thermal_output[h,2:end])
sum_renewable = sum(Renewable_output[h,2:end])
sum_pf = sum_load + sum_thermal + sum_renewable
sum_slackd = sum(slack_dp[h,2:end])
sum_slacku = sum(slack_up[h,2:end])
sum_rload = sum(load_q[h,2:end])
sum_rthermal = sum(Thermal_routput[h,2:end])
sum_rrenewable = sum(Renewable_routput[h,2:end])
sum_slacrkd = sum(slack_dq[h,2:end])
sum_slacrku = sum(slack_uq[h,2:end])

println(sum_rload)
println(sum_rthermal)
println(sum_rrenewable)
println(sum_slacrkd)
println(sum_slacrku)

println(sum_load)
println(sum_thermal)
println(sum_renewable)
println(sum_slackd)
println(sum_slacku)




# Check optimization model
uc1=models.decision_models[1].internal.container
uc2=models.decision_models[2].internal.container

for (k,v) in uc1.constraints
    println(k)
end    
for (k,v) in uc2.constraints
    println(k)
end    

keys(uc1.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{NetworkFlowConstraint, Line}("")] )
uc1.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{NetworkFlowConstraint, Line}("")]["B18",24]


for (k,v) in uc1.variables
    println(k)
end  

for (k,v) in uc2.variables
    println(k)
end  


v=uc2.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowReactivePowerToFromVariable, Line}("")][:,1]  
v=uc2.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowReactivePowerToFromVariable, Line}("")]["B3",1]  

for k in JuMP.all_constraints(uc2.JuMPmodel,; include_variable_in_set_constraints = false)
    #println(k)
    try
        if JuMP.normalized_coefficient(k, v)!=0
            println("coef;",k,";",JuMP.normalized_coefficient(k, v))#, ";is_fixed;",is_fixed(k),";value;",value(k))
            println()
        end
    catch e  end    
end


open("mod_uc1.txt","w") do io
    redirect_stdout(io) do
        println(objective_function(uc1.JuMPmodel))
        for k in all_constraints(uc1.JuMPmodel,; include_variable_in_set_constraints = true)           
            println(name(k),",",k) 
        end    
    end
end
open("mod_uc2.txt","w") do io
    redirect_stdout(io) do
        println(objective_function(uc2.JuMPmodel))
        for k in all_constraints(uc2.JuMPmodel,; include_variable_in_set_constraints = true)           
            println(name(k),",",k) 
        end    
    end
end






