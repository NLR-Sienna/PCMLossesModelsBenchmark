##### Cost-minimization problem: the second stage of the proposed strategy
# This is the final version used to generate the reported results:
# paper_results_version_final, paper_results_version_3_final, and paper_results_version_2_final.
# The solution procedure consists of three steps:
# (1) UC formulated as DCOPF with fixed HVDC power flows (HVDC_value) and fixed losses (Ploss_ori);
# (2) ED formulated as DCOPF with HVDC power flows, a linear loss model, and FND;
# (3) ED formulated as DCOPF with fixed binary variables, HVDC power flows, a linear loss model, and FND.
# FND is incorporated into the transmission constraints in the second and third steps..
# The additional ED step is required because fixing the binary variables obtained
# from the UC model may result in nonzero slack variables in the subsequent ED model.
# Therefore, the binary variables are not fixed in the first ED model.
# Because this model still contains binary variables, valid LMPs cannot be obtained
# directly from the first ED solution.
# Fixed losses are considered in the UC model, while linear losses are considered
# in the ED models.
# The loss factors, estimated power losses, and binary-variable values are obtained
# from the UC model.
# The final outputs include generator dispatch and LMPs.
# Ploss_ori and HVDC_value should be updated for each case.
# Both values are obtained from the first step of the proposed strategy.

using Pkg
Pkg.activate(@__DIR__)
using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using PowerFlows
using TimeSeries
using Ipopt
using Gurobi
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
using JLD2

include("C:/Users/CLI66/.julia/mapped_indices.jl")
include("C:/Users/CLI66/.julia/circular_flows.jl")

const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
const PSI = PowerSimulations
const PNM = PowerNetworkMatrices


##------------------Input data of system---------------------------##
c_sys5_pjm_da = PSB.build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization=true)
T = 1
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(T), Hour(T))

# change the operation cost of renewables
renewables = collect(get_components(RenewableDispatch, c_sys5_pjm_da))
Num_ren = length(renewables)

for row in 1:Num_ren
    new_re_cost = RenewableGenerationCost(;
        variable=CostCurve(LinearCurve(0.01 * row)),
    )
    set_operation_cost!(get_component(RenewableDispatch, c_sys5_pjm_da, get_name(renewables[row])), new_re_cost)
end

hvdc1 = only(get_components(TwoTerminalGenericHVDCLine, c_sys5_pjm_da))
set_loss!(hvdc1, LinearCurve(0.0))

set_reactive_power_limits_from!(hvdc1, (min=-10, max=10))
set_reactive_power_limits_to!(hvdc1, (min=-10, max=10))
set_active_power_limits_from!(hvdc1, (min=-10, max=10))
set_active_power_limits_to!(hvdc1, (min=-10, max=10))

template_uc = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=true, duals=[CopperPlateBalanceConstraint], power_flow_evaluation=PowerFlows.ACPowerFlow(; calculate_loss_factors=true)))
set_device_model!(template_uc, Line, StaticBranchBounds)
set_device_model!(template_uc, TapTransformer, StaticBranchBounds)
set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(
    template_uc,
    DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalDispatch),  # HVDCTwoTerminalLossless
)
set_available!(get_component(HydroDispatch, c_sys5_pjm_da, "201_HYDRO_4"), false)
N_ren=collect(get_components(RenewableNonDispatch, c_sys5_pjm_da))
n_N_ren = length(N_ren)
for j in 1:n_N_ren
    set_available!(N_ren[j], false)
end

#set_rating!(get_component(Line, c_sys5_pjm_da, "C6"), 1.0)
#set_rating!(get_component(Line, c_sys5_pjm_da, "A29"), 0.7)
#set_rating!(get_component(Line, c_sys5_pjm_da, "C33-1"), 0.5)

milp_optimizer = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.001)
nlp_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "MIPGap" => 0.001,
    "OutputFlag" => 1,
    "DisplayInterval" => 1,   # print every node
    "LogFile" => "gurobi.log",
    "Presolve" => 2,
    "Heuristics" => 0.3,
    "NonConvex" => 2,
    "MIPFocus" => 2,        # bound improvement
    "Cuts" => 2
)

UC_models = DecisionModel(
    template_uc,
    c_sys5_pjm_da,
    optimizer=milp_optimizer,
    optimizer_solve_log_print=true,
    name="UC",
    store_variable_names=true,
)

build!(UC_models, output_dir=mktempdir())

# include fixed estimated losses
buses = collect(get_components(ACBus, c_sys5_pjm_da))
num_bus = length(buses)
ref_indexx = PNM.find_slack_positions(buses)
ref_bus = get_number(buses[first(ref_indexx)])
con_bal_uc = UC_models.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][ref_bus, :]
N_con = length(con_bal_uc)
Ploss_ori = 42.2303815088976
HVDC_value = -216.194805651452
base_power = get_base_power(c_sys5_pjm_da)

for j in 1:N_con
    rhs = normalized_rhs(con_bal_uc[j])
    set_normalized_rhs(con_bal_uc[j], rhs + Ploss_ori[j, 1] / base_power)
end

# fix HVDC power flow
HVDC_flow_UC = UC_models.internal.container.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowActivePowerFromToVariable,TwoTerminalGenericHVDCLine}("")]
@constraint(UC_models.internal.container.JuMPmodel, HVDC_flow_UC["DC1", 1] == HVDC_value / 100)

# solve UC and get loss factor, estimated losses and values of binary variables
solve!(UC_models)
optimize!(UC_models.internal.container.JuMPmodel)
results_uc = OptimizationProblemResults(UC_models)
aux_variables_uc = read_aux_variables(results_uc)
# estimated losses
FromTo_Line = aux_variables_uc["PowerFlowLineActivePowerFromTo__Line"]
ToFrom_Line = aux_variables_uc["PowerFlowLineActivePowerToFrom__Line"]
FromTo_TapTransformer = aux_variables_uc["PowerFlowLineActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer = aux_variables_uc["PowerFlowLineActivePowerToFrom__TapTransformer"]
Ploss_est = (sum(Matrix(FromTo_Line)[1, 2:end]) + sum(Matrix(ToFrom_Line)[1, 2:end]) +
    sum(Matrix(FromTo_TapTransformer)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer)[1, 2:end]))
# loss factor
DF = aux_variables_uc["PowerFlowLossFactors__ACBus"] # penalty factor
LF = ones(num_bus, 1) + Matrix(DF)[1, 2:end]
indexx_ref = findfirst(x -> x == ref_bus, parse.(Int, names(DF)[2:end]))
LF[indexx_ref] = 0
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
# power injection
injection_uc = value.(UC_models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
th_power_uc = read_variable(results_uc, "ActivePowerVariable__ThermalStandard")
pcm_uc = objective_value(UC_models.internal.container.JuMPmodel)
re_power_uc = read_variable(results_uc, "ActivePowerVariable__RenewableDispatch")
slack_up_uc = read_variable(results_uc, "SystemBalanceSlackUp__System")
slack_down_uc = read_variable(results_uc, "SystemBalanceSlackDown__System")


## get values of binary variables
ED_models_ini = DecisionModel(
    template_uc,
    c_sys5_pjm_da,
    optimizer=milp_optimizer,
    optimizer_solve_log_print=true,
    name="ED",
    store_variable_names=true,
)

build!(ED_models_ini, output_dir=mktempdir())

# include linear loss
con_bal_ed_ini = ED_models_ini.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][ref_bus, :]
N_con = length(con_bal_ed_ini)
@variable(ED_models_ini.internal.container.JuMPmodel, loss_ed_ini[i=1:N_con, 1], base_name = "loss_ed_ini")
for j in 1:N_con
    set_normalized_coefficient(con_bal_ed_ini[j], loss_ed_ini[j, 1], 1)
    rhs = normalized_rhs(con_bal_ed_ini[j])
    set_normalized_rhs(con_bal_ed_ini[j], rhs + Ploss_est[j, 1] / base_power)
end
injection_ed_ini = ED_models_ini.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")]
AX = axes(injection_ed_ini, 1)
for j in 1:N_con
    @constraint(ED_models_ini.internal.container.JuMPmodel, loss_ed_ini[j, 1] == sum((injection_uc[AX[i], j] - injection_ed_ini[AX[i], j]) * LF[i, j] for i in 1:num_bus))
end


# modify the transmission constraints
con_tran_ed_ini = ED_models_ini.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,Line}("")]
ptdf = PTDF(c_sys5_pjm_da)
for idx in keys(con_tran_ed_ini)
    rhs = normalized_rhs(con_tran_ed_ini[idx])
    indexx = findfirst(x -> x == idx[1], axes(ptdf, 2))
    set_normalized_rhs(con_tran_ed_ini[idx], rhs + sum(ptdf[indexx, i] * E[i, idx[2]] / base_power for i in 1:num_bus))
end

# fix HVDC power flow
HVDC_flow_ini = ED_models_ini.internal.container.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowActivePowerFromToVariable,TwoTerminalGenericHVDCLine}("")]
@constraint(ED_models_ini.internal.container.JuMPmodel, HVDC_flow_ini["DC1", 1] == HVDC_value / 100)

# solve ED and get final generator dispatch and LMP
solve!(ED_models_ini)

optimize!(ED_models_ini.internal.container.JuMPmodel)
results_ed_ini = OptimizationProblemResults(ED_models_ini)
aux_variables_ed_ini = read_aux_variables(results_ed_ini)
# estimated losses
FromTo_Line_ed_ini = aux_variables_ed_ini["PowerFlowLineActivePowerFromTo__Line"]
ToFrom_Line_ed_ini = aux_variables_ed_ini["PowerFlowLineActivePowerToFrom__Line"]
FromTo_TapTransformer_ed_ini = aux_variables_ed_ini["PowerFlowLineActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer_ed_ini = aux_variables_ed_ini["PowerFlowLineActivePowerToFrom__TapTransformer"]
Ploss_est_ed_ini = (sum(Matrix(FromTo_Line_ed_ini)[1, 2:end]) + sum(Matrix(ToFrom_Line_ed_ini)[1, 2:end]) +
    sum(Matrix(FromTo_TapTransformer_ed_ini)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer_ed_ini)[1, 2:end]))
# loss factor
DF_ed_ini = aux_variables_ed_ini["PowerFlowLossFactors__ACBus"] # penalty factor
LF_ed_ini = ones(num_bus, 1) + Matrix(DF_ed_ini)[1, 2:end]
indexx_ref = findfirst(x -> x == ref_bus, parse.(Int, names(DF_ed_ini)[2:end]))
LF_ed_ini[indexx_ref] = 0
# get FND
lines = collect(get_components(Line, c_sys5_pjm_da))
num_line = length(lines)
transformer = collect(get_components(TapTransformer, c_sys5_pjm_da))
num_transformer = length(transformer)
line_loss_ed_ini = zeros(num_line + num_transformer, 1)
for k in 1:num_line
    line_loss_ed_ini[k, 1] = FromTo_Line_ed_ini[1, k+1] + ToFrom_Line_ed_ini[1, k+1]
end
for k in 1:num_transformer
    line_loss_ed_ini[k+num_line, 1] = FromTo_TapTransformer_ed_ini[1, k+1] + ToFrom_TapTransformer_ed_ini[1, k+1]
end
E_ed_ini = zeros(num_bus, 1)
for j in 1:num_line
    arc = get_arc(get_component(Line, c_sys5_pjm_da, get_name(lines[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, names(DF_ed_ini)[2:end]))
    index = findfirst(x -> x == get_name(lines[j]), names(ToFrom_Line_ed_ini)[2:end])
    E_ed_ini[index1, 1] = E_ed_ini[index1, 1] + line_loss_ed_ini[index, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, names(DF_ed_ini)[2:end]))
    E_ed_ini[index2, 1] = E_ed_ini[index2, 1] + line_loss_ed_ini[index, 1] / 2
end
for j in 1:num_transformer
    arc = get_arc(get_component(TapTransformer, c_sys5_pjm_da, get_name(transformer[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, names(DF_ed_ini)[2:end]))
    index = findfirst(x -> x == get_name(transformer[j]), names(ToFrom_TapTransformer_ed_ini)[2:end])
    E_ed_ini[index1, 1] = E_ed_ini[index1, 1] + line_loss_ed_ini[index+num_line, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, names(DF_ed_ini)[2:end]))
    E_ed_ini[index2, 1] = E_ed_ini[index2, 1] + line_loss_ed_ini[index+num_line, 1] / 2
end
# power injection
injection_ed_ini = value.(ED_models_ini.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
th_power_ini = read_variable(results_ed_ini, "ActivePowerVariable__ThermalStandard")
pcm_ini = objective_value(ED_models_ini.internal.container.JuMPmodel)
re_power_ini = read_variable(results_ed_ini, "ActivePowerVariable__RenewableDispatch")
slack_up_ini = read_variable(results_ed_ini, "SystemBalanceSlackUp__System")
slack_down_ini = read_variable(results_ed_ini, "SystemBalanceSlackDown__System")


## run ED with linear loss and fixed HVDC power flow
ED_models = DecisionModel(
    template_uc,
    c_sys5_pjm_da,
    optimizer=milp_optimizer,
    optimizer_solve_log_print=true,
    name="ED",
    store_variable_names=true,
)

build!(ED_models, output_dir=mktempdir())

# fix binary variables based on the results from UC_models
fixed_uc_variables = all_variables(ED_models_ini.internal.container.JuMPmodel)
fixed_uc_solution_by_name = Dict(zip(name.(fixed_uc_variables), value.(fixed_uc_variables)))
ed_variables = all_variables(ED_models.internal.container.JuMPmodel)
for ed_var in ed_variables
    if is_binary(ed_var)
        unset_binary(ed_var)
        set_upper_bound(ed_var, fixed_uc_solution_by_name[name(ed_var)])
        set_lower_bound(ed_var, fixed_uc_solution_by_name[name(ed_var)])
    end
end

# include linear loss
con_bal_ed = ED_models.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][ref_bus, :]
N_con = length(con_bal_ed)
@variable(ED_models.internal.container.JuMPmodel, loss_ed[i=1:N_con, 1], base_name = "loss_ed")
for j in 1:N_con
    set_normalized_coefficient(con_bal_ed[j], loss_ed[j, 1], 1)
    rhs = normalized_rhs(con_bal_ed[j])
    set_normalized_rhs(con_bal_ed[j], rhs + Ploss_est_ed_ini[j, 1] / base_power)
end
injection_ed = ED_models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")]
AX = axes(injection_ed, 1)
for j in 1:N_con
    @constraint(ED_models.internal.container.JuMPmodel, loss_ed[j, 1] == sum((injection_ed_ini[AX[i], j] - injection_ed[AX[i], j]) * LF_ed_ini[i, j] for i in 1:num_bus))
end

# modify the transmission constraints
con_tran_ed = ED_models.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,Line}("")]
ptdf = PTDF(c_sys5_pjm_da)
for idx in keys(con_tran_ed)
    rhs = normalized_rhs(con_tran_ed[idx])
    indexx = findfirst(x -> x == idx[1], axes(ptdf, 2))
    set_normalized_rhs(con_tran_ed[idx], rhs + sum(ptdf[indexx, i] * E_ed_ini[i, idx[2]] / base_power for i in 1:num_bus))
end

# fix HVDC power flow
HVDC_flow = ED_models.internal.container.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowActivePowerFromToVariable,TwoTerminalGenericHVDCLine}("")]
@constraint(ED_models.internal.container.JuMPmodel, HVDC_flow["DC1", 1] == HVDC_value / 100)

# solve ED and get final generator dispatch and LMP
solve!(ED_models)
optimize!(ED_models.internal.container.JuMPmodel)
pcm = objective_value(ED_models.internal.container.JuMPmodel)
results_ed = OptimizationProblemResults(ED_models)
th_power_ed = read_variable(results_ed, "ActivePowerVariable__ThermalStandard")
load_param_ed = read_parameter(results_ed, "ActivePowerTimeSeriesParameter__PowerLoad")
renewable_limt = read_parameter(results_ed, "ActivePowerTimeSeriesParameter__RenewableDispatch")
re_power_ed = read_variable(results_ed, "ActivePowerVariable__RenewableDispatch")
slack_up_ed = read_variable(results_ed, "SystemBalanceSlackUp__System")
slack_down_ed = read_variable(results_ed, "SystemBalanceSlackDown__System")
injection_1 = value.(ED_models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
loss_0 = zeros(T, 1)
for j in 1:T
    loss_0[j, 1] = sum(injection_1[:, j])
end

# get power losses based on AC power flow calculation
aux_variables_ed = read_aux_variables(results_ed)
FromTo_Line_ed = aux_variables_ed["PowerFlowLineActivePowerFromTo__Line"]
ToFrom_Line_ed = aux_variables_ed["PowerFlowLineActivePowerToFrom__Line"]
FromTo_TapTransformer_ed = aux_variables_ed["PowerFlowLineActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer_ed = aux_variables_ed["PowerFlowLineActivePowerToFrom__TapTransformer"]
Ploss_ed = (sum(Matrix(FromTo_Line_ed)[1, 2:end]) + sum(Matrix(ToFrom_Line_ed)[1, 2:end]) +
    sum(Matrix(FromTo_TapTransformer_ed)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer_ed)[1, 2:end]))


# get line flow from DCOPF calculation
flow_line_1 = read_variable(results_ed, "FlowActivePowerVariable__Line")
flow_transformer_1 = read_variable(results_ed, "FlowActivePowerVariable__TapTransformer")

HVDC_flow_1 = read_variable(results_ed, "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine")
HVDC_loss = read_variable(results_ed, "HVDCLosses__TwoTerminalGenericHVDCLine")


# calculate LMP
lamda = read_dual(results_ed, "CopperPlateBalanceConstraint__System") # energy price
miu = zeros(1, length(axes(ptdf, 2)))
net_con = ED_models.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,Line}("")]
net_con_name = axes(net_con, 1)
net_con2 = ED_models.internal.container.constraints[PowerSimulations.ConstraintKey{NetworkFlowConstraint,TapTransformer}("")]
net_con_name2 = axes(net_con2, 1)
for j in 1:(num_line+num_transformer)
    index1 = findfirst(x -> x == axes(ptdf, 2)[j], net_con_name)
    if isnothing(index1)
        index2 = findfirst(x -> x == axes(ptdf, 2)[j], net_con_name2)
        miu[1, j] = dual(net_con2[net_con_name2[index2, 1], 1])
    else
        miu[1, j] = dual(net_con[net_con_name[index1, 1], 1])
    end
end
LMP = zeros(num_bus, 1)
for j in 1:num_bus
    LMP[j, 1] = lamda[1, 2] + (miu*ptdf[:, j])[1, 1] + lamda[1, 2] * (-LF_ed_ini[j, 1])
end

# calculate load payment
load_bus_index = names(load_param_ed)[2:end]
load_bus_num = length(load_bus_index)
load_payment = zeros(load_bus_num, 1)
for j in 1:load_bus_num
    indexx = findfirst(x -> x == get_number(get_component(ACBus, c_sys5_pjm_da, load_bus_index[j])), axes(ptdf, 1))
    load_payment[j, 1] = -load_param_ed[1, j+1]*LMP[indexx, 1]
end
sum_load_payment = sum(load_payment)

# calculate generator revenue
thermal_generator = collect(get_components(ThermalStandard, c_sys5_pjm_da))
num_thermal = length(thermal_generator)
thermal_revenue = zeros(num_thermal, 1)
renewable_revenue = zeros(Num_ren, 1)
names_thermal = names(th_power_ed)[2:end]
names_renewable = names(re_power_ed)[2:end]
for j in 1:num_thermal
    indexx = findfirst(x -> x == get_number(get_bus(get_component(ThermalStandard, c_sys5_pjm_da, names_thermal[j]))), axes(ptdf, 1))
    thermal_revenue[j, 1] = th_power_ed[1, j+1]*LMP[indexx, 1]
end
for j in 1:Num_ren
    indexx = findfirst(x -> x == get_number(get_bus(get_component(RenewableDispatch, c_sys5_pjm_da, names_renewable[j]))), axes(ptdf, 1))
    renewable_revenue[j, 1] = re_power_ed[1, j+1]*LMP[indexx, 1]
end
sum_generator_revenue = sum(thermal_revenue) + sum(renewable_revenue)
renewable_subsidy = 0.3 * 100 * sum(Matrix(re_power_ed)[2:end])

(sum_load_payment-sum_generator_revenue)/lamda[1, 2]
sum(injection_ed_ini[AX[i], 1] * LF_ed_ini[i, 1] for i in 1:num_bus)*100 - Ploss_est_ed_ini



#=
# output the final results
XLSX.openxlsx("Results_KKT_2_A29P.xlsx", mode="w") do xf
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "results")
    sheet1["A1"] = "obj"
    sheet1["A2"] = pcm
    sheet1["B1"] = "AC_loss"
    sheet1["B2"] = Ploss_ed
    sheet1["C1"] = "DC_loss"
    sheet1["C2"] = loss_0[1]
    sheet1["D1"] = "calculate_loss"
    sheet1["D2"] = value(loss_ed[1, 1])
    sheet1["E1"] = "HVDC_flow"
    sheet1["E2"] = HVDC_flow_1[1, 2]

    sheet2 = XLSX.addsheet!(xf, "Generator")
    XLSX.writetable!(sheet2, Tables.table(Matrix(th_power_ed)[2:end]))
    sheet3 = XLSX.addsheet!(xf, "Renewable_Gen")
    XLSX.writetable!(sheet3, Tables.table(Matrix(re_power_ed)[2:end]))
    sheet4 = XLSX.addsheet!(xf, "LMP")
    XLSX.writetable!(sheet4, Tables.table(LMP))
    sheet5 = XLSX.addsheet!(xf, "line_flow")
    XLSX.writetable!(sheet5, Tables.table(Matrix(flow_line_1)[2:end]))
    sheet6 = XLSX.addsheet!(xf, "transformer_flow")
    XLSX.writetable!(sheet6, Tables.table(Matrix(flow_transformer_1)[2:end]))
end
=#

# identify cycles
data = PSI.get_power_flow_data(
    only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(ED_models))),
)
base_power = get_base_power(c_sys5_pjm_da)

res_vars = read_variables(results_ed)
time_step = 1
G = build_graph(data; time_step=time_step)
add_hvdc_edges!(G, c_sys5_pjm_da, res_vars, data; time_step=time_step)
branches = collect(get_components(ACBranch, c_sys5_pjm_da))
C = find_circular_flows(G, data, branches)
num_cycle = length(C)
circular_flow = []
for j in 1:num_cycle
    if "DC1" in C[j].branches
        push!(circular_flow, minimum(C[j, 1].branch_flows))
    end
end


XLSX.openxlsx("circular flow.xlsx", mode="w") do xf
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "cycles")
    XLSX.writetable!(sheet1, Tables.table(circular_flow))
end

println(slack_up_ed)
println(slack_down_ed)
println(C)


#=
uc0 = ED_models.internal.container
objd0 = Dict()
sol0 = Dict()
objv0 = 0
obj0 = objective_function(uc0.JuMPmodel)
objd0 = Dict(zip(name.(keys(obj0.terms)), values(obj0.terms)))
av0 = JuMP.all_variables(uc0.JuMPmodel)
sol0 = Dict(zip(name.(av0), value.(av0)))
namess = []
valuess = []
objcoef = []
obj_contribution = []
for (k, v) in objd0
    objv0 = objv0 + sol0[k] * objd0[k]
    #println("name;value;objcoef;obj_contribution;", k, ";", sol0[k], ";", objd0[k], ";", sol0[k] * objd0[k])
    push!(namess, k)
    push!(valuess, sol0[k])
    push!(objcoef, objd0[k])
    push!(obj_contribution, sol0[k] * objd0[k])
end

XLSX.openxlsx("objective_term.xlsx", mode="w") do xf
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "name")
    XLSX.writetable!(sheet1, Tables.table(namess))
    sheet2 = XLSX.addsheet!(xf, "value")
    XLSX.writetable!(sheet2, Tables.table(valuess))
    sheet3 = XLSX.addsheet!(xf, "objcoef")
    XLSX.writetable!(sheet3, Tables.table(objcoef))
    sheet4 = XLSX.addsheet!(xf, "obj_contribution")
    XLSX.writetable!(sheet4, Tables.table(obj_contribution))
end
=#