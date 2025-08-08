# solve the SCUC problem considering linear loss for RTS
# Step 1: solve DCOPF without losses
# Step 2: perform iteration based on AC power flow and FND
# consider renewable generators
# do not consider HVDC
# output LMP
# power balance constraint is: (G-D) + LF*(injection2-injection1)=P_loss^est

using Pkg
Pkg.activate(@__DIR__)
using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using PowerFlows
using TimeSeries
#using Ipopt
using HiGHS
using Logging
using Revise
using JuMP
using Dates
using CSV
using DataFrames
using HydroPowerSimulations
using PowerNetworkMatrices
using XLSX


const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
const PSI = PowerSimulations
const PNM = PowerNetworkMatrices


##------------------Input data of system---------------------------##
c_sys5_pjm_da = PSB.build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization=true)
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(1), Hour(1))

template_uc = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=true, duals=[CopperPlateBalanceConstraint], power_flow_evaluation=PowerFlows.ACPowerFlow(; calculate_loss_factors=true)))
set_device_model!(template_uc, Line, StaticBranchBounds)
set_device_model!(template_uc, TapTransformer, StaticBranchBounds)
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, TwoTerminalHVDCLine, HVDCTwoTerminalDispatch)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)

milp_optimizer = optimizer_with_attributes(HiGHS.Optimizer)

set_available!(get_component(TwoTerminalHVDCLine, c_sys5_pjm_da, "DC1"), false)
#set_rating!(get_component(Line, c_sys5_pjm_da, "A23"), 2.5)

models = DecisionModel(
    template_uc,
    c_sys5_pjm_da,
    optimizer=milp_optimizer,
    name="UC",
    store_variable_names=true,
)

build!(models, output_dir=mktempdir())

##------------------Step 0,1 solve UC without losses and output results---------------------------##
solve!(models)

optimize!(models.internal.container.JuMPmodel)
pcm_1 = objective_value(models.internal.container.JuMPmodel)
# obtain the injection at each bus
injection_1 = value.(models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])

## Output results
# Retriving results
results_uc_1 = OptimizationProblemResults(models)
th_power_1 = read_variable(results_uc_1, "ActivePowerVariable__ThermalStandard")
load_param_1 = read_parameter(results_uc_1, "ActivePowerTimeSeriesParameter__PowerLoad")
renewable_limt = read_parameter(results_uc_1, "ActivePowerTimeSeriesParameter__RenewableDispatch")
re_power_1 = read_variable(results_uc_1, "ActivePowerVariable__RenewableDispatch")
hy_power_1 = read_variable(results_uc_1, "ActivePowerVariable__HydroDispatch")

slack_up_1 = read_variable(results_uc_1, "SystemBalanceSlackUp__System")
slack_down_1 = read_variable(results_uc_1, "SystemBalanceSlackDown__System")
loss_0 = sum(Matrix(th_power_1)[1, 2:end]) + sum(Matrix(re_power_1)[1, 2:end]) + hy_power_1[1, 2] + sum(Matrix(load_param_1)[1, 2:end])

aux_variables = read_aux_variables(results_uc_1)
# get loss factor based on AC power flow calculation
DF = aux_variables["PowerFlowLossFactors__ACBus"] # penalty factor
buses = collect(get_components(ACBus, c_sys5_pjm_da))
num_bus = length(buses)
LF = ones(num_bus, 1) + Matrix(DF)[1, 2:end]
indexx = PNM.find_slack_positions(buses)
ref_bus = get_number(buses[first(indexx)])
index = findfirst(x -> x == ref_bus, parse.(Int, names(DF)[2:end]))
LF[index] = 0

# get estimated power losses based on AC power flow calculation
FromTo_Line = aux_variables["PowerFlowLineActivePowerFromTo__Line"]
ToFrom_Line = aux_variables["PowerFlowLineActivePowerToFrom__Line"]
FromTo_TapTransformer = aux_variables["PowerFlowLineActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer = aux_variables["PowerFlowLineActivePowerToFrom__TapTransformer"]
Ploss_est = (sum(Matrix(FromTo_Line)[1, 2:end]) + sum(Matrix(ToFrom_Line)[1, 2:end]) +
             sum(Matrix(FromTo_TapTransformer)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer)[1, 2:end]))

# get FND
lines = collect(get_components(Line, c_sys5_pjm_da))
num_line = length(lines)
transformer = collect(get_components(TapTransformer, c_sys5_pjm_da))
num_transformer = length(transformer)
line_loss = zeros(num_line + num_transformer, 1)
for k in 1:num_line
    line_loss[k, 1] = FromTo_Line[1, k+1] + ToFrom_Line[1, k+1]
end
for k in 1:num_transformer
    line_loss[k+num_line, 1] = FromTo_TapTransformer[1, k+1] + ToFrom_TapTransformer[1, k+1]
end

E = zeros(num_bus, 1)
for j in 1:num_line
    arc = get_arc(get_component(Line, c_sys5_pjm_da, get_name(lines[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, names(DF)[2:end]))
    index = findfirst(x -> x == get_name(lines[j]), names(ToFrom_Line)[2:end])
    E[index1, 1] = E[index1, 1] + line_loss[index, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, names(DF)[2:end]))
    E[index2, 1] = E[index2, 1] + line_loss[index, 1] / 2
end
for j in 1:num_transformer
    arc = get_arc(get_component(TapTransformer, c_sys5_pjm_da, get_name(transformer[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, names(DF)[2:end]))
    index = findfirst(x -> x == get_name(transformer[j]), names(ToFrom_TapTransformer)[2:end])
    E[index1, 1] = E[index1, 1] + line_loss[index+num_line, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, names(DF)[2:end]))
    E[index2, 1] = E[index2, 1] + line_loss[index+num_line, 1] / 2
end

# get line flow from DCOPF calculation
flow_line_1 = read_variable(results_uc_1, "FlowActivePowerVariable__Line")
flow_transformer_1 = read_variable(results_uc_1, "FlowActivePowerVariable__TapTransformer")

# get LMP, lamda is the energy price, miu is the dual of transmission constraints
ptdf = PTDF(c_sys5_pjm_da)
lamda = read_dual(results_uc_1, "CopperPlateBalanceConstraint__System") # energy price
miu = zeros(1, length(axes(ptdf, 2)))
net_con = models.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,Line}("")]
net_con_name = axes(net_con, 1)
net_con2 = models.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,TapTransformer}("")]
net_con_name2 = axes(net_con2, 1)
for j in 1:num_line+num_transformer
    index1 = findfirst(x -> x == axes(ptdf, 2)[j], net_con_name)
    if isnothing(index1)
        index2 = findfirst(x -> x == axes(ptdf, 2)[j], net_con_name2)
        miu[1, j] = dual(net_con2[net_con_name2[index2, 1], 1])
    else
        miu[1, j] = dual(net_con[net_con_name[index1, 1], 1])
    end
end
LMP_1 = zeros(num_bus, 1)
for j in 1:num_bus
    LMP_1[j, 1] = lamda[1, 2] + (miu*ptdf[:, j])[1, 1]
end

##------------------Step 2 perfomr iteration to consider linear loss in SCUC problem---------------------------##
iter_max = 10 # set the maximum number of iteration
# set the initial value of the result to be saved
num_gen = length(Matrix(th_power_1)) - 1
th_res = zeros(iter_max + 1, num_gen)
re_res = zeros(iter_max + 1, Num_ren)
LMP = zeros(num_bus, iter_max + 1)
LMP[:, 1] = LMP_1
err = zeros(iter_max + 1, 1)
th_res[1, :] = Matrix(th_power_1)[1, 2:1+num_gen]
re_res[1, :] = Matrix(re_power_1)[1, 2:1+Num_ren]
err[1, 1] = 10
err_sta = 0.01
base_power = get_base_power(c_sys5_pjm_da)
ptdf = PTDF(c_sys5_pjm_da)
LOSS_DC = zeros(iter_max + 1, 1)#DC loss obtained from DCOPF
PLOSS = zeros(iter_max + 1, 1) #AC loss obtained from AC power flow calculation
PLOSS[1, 1] = Ploss_est
PCM = zeros(iter_max + 1, 1)
PCM[1, 1] = pcm_1
LF_iter = zeros(num_bus, iter_max + 1)
LF_iter[:, 1] = LF
Injec = zeros(num_bus, iter_max + 1)
Injec[:, 1] = Matrix(injection_1)

line_flow_AC = zeros(num_line, iter_max + 1)
line_flow_DC = zeros(num_line, iter_max + 1)
transformer_flow_AC = zeros(num_transformer, iter_max + 1)
transformer_flow_DC = zeros(num_transformer, iter_max + 1)
line_flow_AC[:, 1] = Matrix(FromTo_Line)[1, 2:end]
transformer_flow_AC[:, 1] = Matrix(FromTo_TapTransformer)[1, 2:end]
line_flow_DC[:, 1] = Matrix(flow_line_1)[1, 2:end]
transformer_flow_DC[:, 1] = Matrix(flow_transformer_1)[1, 2:end]

for iter in 1:iter_max
    global ptdf
    global Ploss_est
    global LF
    global E
    global injection_1

    models_2 = DecisionModel(
        template_uc,
        c_sys5_pjm_da,
        optimizer=milp_optimizer,
        name="UC_2",
        store_variable_names=true,
    )

    build!(models_2, output_dir=mktempdir())

    # add losses in energy balance constraint
    NT = length(models_2.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][ref_bus, :])
    con_bal = models_2.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][ref_bus, :]
    @variable(models_2.internal.container.JuMPmodel, llosses[i=1:NT, 1], base_name = "line_losses")

    # include losses in the energy balance constraint
    num_con_bal = length(con_bal)
    for j in 1:num_con_bal
        set_normalized_coefficient(con_bal[j], llosses[j, 1], 1)
        rhs = normalized_rhs(con_bal[j])
        set_normalized_rhs(con_bal[j], rhs + Ploss_est[j, 1] / base_power)
    end

    # calculate the power losses based on loss factor
    injection = models_2.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")]
    AX = axes(injection, 1)

    for j in 1:NT
        @constraint(models_2.internal.container.JuMPmodel, llosses[j, 1] == sum((injection_1[AX[i], j] - injection[AX[i], j]) * LF[i, j] for i in 1:num_bus))
    end

    ## modify the transmission constraints
    con_tran = models_2.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,Line}("")]
    ptdf = PTDF(c_sys5_pjm_da)
    for idx in keys(con_tran)
        rhs = normalized_rhs(con_tran[idx])
        indexx = findfirst(x -> x == idx[1], axes(ptdf, 2))
        set_normalized_rhs(con_tran[idx], rhs + sum(ptdf[indexx, i] * E[i, idx[2]] / base_power for i in 1:num_bus))
    end
    solve!(models_2)

    ## Output results
    # Retriving results
    results_uc_2 = OptimizationProblemResults(models_2)

    optimize!(models_2.internal.container.JuMPmodel)
    pcm_2 = objective_value(models_2.internal.container.JuMPmodel)
    load_param_2 = read_parameter(results_uc_2, "ActivePowerTimeSeriesParameter__PowerLoad")
    th_power_2 = read_variable(results_uc_2, "ActivePowerVariable__ThermalStandard")
    re_power_2 = read_variable(results_uc_2, "ActivePowerVariable__RenewableDispatch")
    hy_power_2 = read_variable(results_uc_2, "ActivePowerVariable__HydroDispatch")
    slack_up_2 = read_variable(results_uc_2, "SystemBalanceSlackUp__System")
    slack_down_2 = read_variable(results_uc_2, "SystemBalanceSlackDown__System")
    loss_DC_2 = sum(Matrix(th_power_2)[1, 2:end]) + sum(Matrix(re_power_2)[1, 2:end]) + hy_power_2[1, 2] + sum(Matrix(load_param_2)[1, 2:end])

    injection_1 = value.(models_2.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
    injection_2 = value.(models_2.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
    LF_iter[:, iter+1] = LF

    # get LMP
    lamda = read_dual(results_uc_2, "CopperPlateBalanceConstraint__System") # energy price
    miu = zeros(1, length(axes(ptdf, 2))) # dual variable of transmission constraint
    net_con = models_2.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,Line}("")]
    net_con_name = axes(net_con, 1)
    net_con2 = models_2.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,TapTransformer}("")]
    net_con_name2 = axes(net_con2, 1)
    for j in 1:num_line+num_transformer
        index1 = findfirst(x -> x == axes(ptdf, 2)[j], net_con_name)
        if isnothing(index1)
            index2 = findfirst(x -> x == axes(ptdf, 2)[j], net_con_name2)
            miu[1, j] = dual(net_con2[net_con_name2[index2, 1], 1])
        else
            miu[1, j] = dual(net_con[net_con_name[index1, 1], 1])
        end
    end
    LMP_2 = zeros(num_bus, 1)
    for j in 1:num_bus
        LMP_2[j, 1] = lamda[1, 2] + (miu*ptdf[:, j])[1, 1] + lamda[1, 2] * (-LF[j])
    end
    LMP[:, iter+1] = LMP_2

    # get line flow, power loss, and loss factor from AC power flow calculation
    aux_variables_2 = read_aux_variables(results_uc_2)
    FromTo_Line_2 = aux_variables_2["PowerFlowLineActivePowerFromTo__Line"]
    ToFrom_Line_2 = aux_variables_2["PowerFlowLineActivePowerToFrom__Line"]
    FromTo_TapTransformer_2 = aux_variables_2["PowerFlowLineActivePowerFromTo__TapTransformer"]
    ToFrom_TapTransformer_2 = aux_variables_2["PowerFlowLineActivePowerToFrom__TapTransformer"]
    Ploss_est = (sum(Matrix(FromTo_Line_2)[1, 2:end]) + sum(Matrix(ToFrom_Line_2)[1, 2:end]) +
                 sum(Matrix(FromTo_TapTransformer_2)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer_2)[1, 2:end]))
    DF = aux_variables_2["PowerFlowLossFactors__ACBus"] # penalty factor               
    LF = ones(num_bus, 1) + Matrix(DF)[1, 2:end]
    index = findfirst(x -> x == ref_bus, parse.(Int, names(DF)[2:end]))
    LF[index] = 0

    # calculate FND
    line_loss = zeros(num_line + num_transformer, 1)
    for k in 1:num_line
        line_loss[k, 1] = FromTo_Line_2[1, k+1] + ToFrom_Line_2[1, k+1]
    end
    for k in 1:num_transformer
        line_loss[k+num_line, 1] = FromTo_TapTransformer_2[1, k+1] + ToFrom_TapTransformer_2[1, k+1]
    end

    E = zeros(num_bus, 1)
    for j in 1:num_line
        arc = get_arc(get_component(Line, c_sys5_pjm_da, get_name(lines[j])))
        bus1 = get_from(arc)
        num1 = get_number(bus1)
        index1 = findfirst(x -> x == num1, parse.(Int, names(DF)[2:end]))
        index = findfirst(x -> x == get_name(lines[j]), names(ToFrom_Line_2)[2:end])
        E[index1, 1] = E[index1, 1] + line_loss[index, 1] / 2
        bus2 = get_to(arc)
        num2 = get_number(bus2)
        index2 = findfirst(x -> x == num2, parse.(Int, names(DF)[2:end]))
        E[index2, 1] = E[index2, 1] + line_loss[index, 1] / 2
    end
    for j in 1:num_transformer
        arc = get_arc(get_component(TapTransformer, c_sys5_pjm_da, get_name(transformer[j])))
        bus1 = get_from(arc)
        num1 = get_number(bus1)
        index1 = findfirst(x -> x == num1, parse.(Int, names(DF)[2:end]))
        index = findfirst(x -> x == get_name(transformer[j]), names(ToFrom_TapTransformer_2)[2:end])
        E[index1, 1] = E[index1, 1] + line_loss[index+num_line, 1] / 2
        bus2 = get_to(arc)
        num2 = get_number(bus2)
        index2 = findfirst(x -> x == num2, parse.(Int, names(DF)[2:end]))
        E[index2, 1] = E[index2, 1] + line_loss[index+num_line, 1] / 2
    end

    # get the line flow from DCOPF
    flow_line_2 = read_variable(results_uc_2, "FlowActivePowerVariable__Line")
    flow_transformer_2 = read_variable(results_uc_2, "FlowActivePowerVariable__TapTransformer")

    # save results
    line_flow_DC[:, iter+1] = Matrix(flow_line_2)[1, 2:end]
    transformer_flow_DC[:, iter+1] = Matrix(flow_transformer_2)[1, 2:end]
    line_flow_AC[:, iter+1] = Matrix(FromTo_Line_2)[1, 2:end]
    transformer_flow_AC[:, iter+1] = Matrix(FromTo_TapTransformer_2)[1, 2:end]

    Injec[:, iter+1] = Matrix(injection_2)
    PLOSS[iter+1, 1] = Ploss_est
    LOSS_DC[iter+1, 1] = loss_DC_2
    th_res[iter+1, :] = Matrix(th_power_2)[1, 2:1+num_gen]
    re_res[iter+1, :] = Matrix(re_power_2)[1, 2:1+Num_ren]
    PCM[iter+1, 1] = pcm_2

    # use the difference in the generators' output between two iterations as the criterion for judging convergence
    err1 = maximum(abs.(th_res[iter+1, :] - th_res[iter, :]))
    err2 = maximum(abs.(re_res[iter+1, :] - re_res[iter, :]))
    err3 = maximum([err1, err2])
    err[iter+1, 1] = err3

    # use the difference in the objective function between two iterations as the criterion for judging convergence
    #err[iter+1, 1] = abs(PCM[iter+1, 1] - PCM[iter, 1]) / PCM[iter, 1] * 100
    if err[iter+1, 1] <= err_sta
        # output final LMP after convergence
        XLSX.openxlsx("LMP.xlsx", mode="w") do xf
            sheet1 = xf[1]
            #sheet1 = XLSX.addsheet!(xf, "LMP")
            XLSX.rename!(sheet1, "LMP")
            XLSX.writetable!(sheet1, Tables.table(LMP_2))

            sheet2 = XLSX.addsheet!(xf, "energy_price")
            sheet2["A1"] = "price"
            sheet2["A2"] = lamda[1, 2]

            sheet3 = XLSX.addsheet!(xf, "congestion_price")
            XLSX.writetable!(sheet3, Tables.table(miu'))

            sheet4 = XLSX.addsheet!(xf, "line_flow")
            XLSX.writetable!(sheet4, Tables.table(Matrix(flow_line_2)[1, 2:end]))
        end
        break
    end
end

# output the final results
XLSX.openxlsx("Results.xlsx", mode="w") do xf
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "results")
    sheet1["A1"] = "obj"
    sheet1["A2:A12"] = PCM
    sheet1["B1"] = "AC_loss"
    sheet1["B2:B12"] = PLOSS
    sheet1["C1"] = "DC_loss"
    sheet1["C2:C12"] = LOSS_DC

    sheet2 = XLSX.addsheet!(xf, "Thermal_gen")
    XLSX.writetable!(sheet2, Tables.table(th_res))
    sheet3 = XLSX.addsheet!(xf, "Renewable_Gen")
    XLSX.writetable!(sheet3, Tables.table(re_res))
    sheet4 = XLSX.addsheet!(xf, "LMP")
    XLSX.writetable!(sheet4, Tables.table(LMP))
    sheet5 = XLSX.addsheet!(xf, "LF_iter")
    XLSX.writetable!(sheet5, Tables.table(LF_iter))
    sheet6 = XLSX.addsheet!(xf, "line_flow_AC")
    XLSX.writetable!(sheet6, Tables.table(line_flow_AC))
    sheet7 = XLSX.addsheet!(xf, "transformer_flow_AC")
    XLSX.writetable!(sheet7, Tables.table(transformer_flow_AC))
    sheet8 = XLSX.addsheet!(xf, "line_flow_DC")
    XLSX.writetable!(sheet8, Tables.table(line_flow_DC))
    sheet9 = XLSX.addsheet!(xf, "transformer_flow_DC")
    XLSX.writetable!(sheet9, Tables.table(transformer_flow_DC))
    sheet10 = XLSX.addsheet!(xf, "Injec")
    XLSX.writetable!(sheet10, Tables.table(Injec))
end


