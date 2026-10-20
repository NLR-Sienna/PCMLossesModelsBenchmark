
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

include("mapped_indices.jl") ## should include this file
include("circular_flows_V2.jl") ## should include this file

const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
const PSI = PowerSimulations
const PNM = PowerNetworkMatrices

LOCAL_PATH = pwd()
steps_sim = 1
folder_path = joinpath(LOCAL_PATH, "Sienna Results", "$(today())_simple_psse_export_$(steps_sim)steps")
mkpath(folder_path)
psse_export_dir = joinpath(folder_path, "psse_export")
mkpath(psse_export_dir)

##------------------Input data of system---------------------------##
c_sys5_pjm_da = PSB.build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization=true)
T = 1
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(T), Hour(T))

# change the operation cost of renewables
renewables = collect(get_components(RenewableDispatch, c_sys5_pjm_da))
Num_ren = length(renewables)

for row in 1:Num_ren
    new_re_cost = RenewableGenerationCost(;
        variable=CostCurve(LinearCurve(-0.01 * row)),
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

milp_optimizer = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.001)

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
FromTo_Line = aux_variables_uc["PowerFlowBranchActivePowerFromTo__Line"]
ToFrom_Line = aux_variables_uc["PowerFlowBranchActivePowerToFrom__Line"]
FromTo_TapTransformer = aux_variables_uc["PowerFlowBranchActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer = aux_variables_uc["PowerFlowBranchActivePowerToFrom__TapTransformer"]
Ploss_est = (sum(Matrix(FromTo_Line)[:, 3]) + sum(Matrix(ToFrom_Line)[:, 3]) +
    sum(Matrix(FromTo_TapTransformer)[:, 3]) + sum(Matrix(ToFrom_TapTransformer)[:, 3]))

line_loss = aux_variables_uc["PowerFlowBranchActivePowerLoss__Line"]
tap_loss = aux_variables_uc["PowerFlowBranchActivePowerLoss__TapTransformer"]
#HVDC_loss = aux_variables_uc["PowerFlowBranchActivePowerLoss__TwoTerminalGenericHVDCLine"]

# loss factor
DF = aux_variables_uc["PowerFlowLossFactors__ACBus"] # penalty factor
LF = ones(num_bus, 1) + Matrix(DF)[:, 3]
indexx_ref = findfirst(x -> x == ref_bus, parse.(Int, DF[:, 2]))
LF[indexx_ref] = 0
# get FND
lines = collect(get_components(Line, c_sys5_pjm_da))
num_line = length(lines)
transformer = collect(get_components(TapTransformer, c_sys5_pjm_da))
num_transformer = length(transformer)
line_loss = zeros(num_line + num_transformer, 1)
for k in 1:num_line
    line_loss[k, 1] = FromTo_Line[k, 3] + ToFrom_Line[k, 3]
end
for k in 1:num_transformer
    line_loss[k+num_line, 1] = FromTo_TapTransformer[k, 3] + ToFrom_TapTransformer[k, 3]
end
E = zeros(num_bus, 1)
for j in 1:num_line
    arc = get_arc(get_component(Line, c_sys5_pjm_da, get_name(lines[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, DF[:, 2]))
    index = findfirst(x -> x == get_name(lines[j]), ToFrom_Line[:, 2])
    E[index1, 1] = E[index1, 1] + line_loss[index, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, DF[:, 2]))
    E[index2, 1] = E[index2, 1] + line_loss[index, 1] / 2
end
for j in 1:num_transformer
    arc = get_arc(get_component(TapTransformer, c_sys5_pjm_da, get_name(transformer[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, DF[:, 2]))
    index = findfirst(x -> x == get_name(transformer[j]), ToFrom_TapTransformer[:, 2])
    E[index1, 1] = E[index1, 1] + line_loss[index+num_line, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, DF[:, 2]))
    E[index2, 1] = E[index2, 1] + line_loss[index+num_line, 1] / 2
end
# power injection
injection_uc = value.(UC_models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])

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
    line_name = idx[1]
    if endswith(line_name, "-double_circuit")
        line_name = replace(line_name, r"-double_circuit$" => "-1")
    end
    line = get_component(Line, c_sys5_pjm_da, line_name)
    from = get_number(get_from(get_arc(line)))
    to = get_number(get_to(get_arc(line)))
    indexx = findfirst(
        x -> x == (from, to),
        axes(ptdf, 2)
    )
    set_normalized_rhs(
        con_tran_ed_ini[idx],
        rhs + sum(
            ptdf[indexx, i] * E[i, 1] / base_power
            for i in 1:num_bus
        )
    )
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
FromTo_Line_ed_ini = aux_variables_ed_ini["PowerFlowBranchActivePowerFromTo__Line"]
ToFrom_Line_ed_ini = aux_variables_ed_ini["PowerFlowBranchActivePowerToFrom__Line"]
FromTo_TapTransformer_ed_ini = aux_variables_ed_ini["PowerFlowBranchActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer_ed_ini = aux_variables_ed_ini["PowerFlowBranchActivePowerToFrom__TapTransformer"]
Ploss_est_ed_ini = (sum(Matrix(FromTo_Line_ed_ini)[:, 3]) + sum(Matrix(ToFrom_Line_ed_ini)[:, 3]) +
    sum(Matrix(FromTo_TapTransformer_ed_ini)[:, 3]) + sum(Matrix(ToFrom_TapTransformer_ed_ini)[:, 3]))
# loss factor
DF_ed_ini = aux_variables_ed_ini["PowerFlowLossFactors__ACBus"] # penalty factor
LF_ed_ini = ones(num_bus, 1) + Matrix(DF_ed_ini)[:, 3]
indexx_ref = findfirst(x -> x == ref_bus, parse.(Int, DF_ed_ini[:, 2]))
LF_ed_ini[indexx_ref] = 0
# get FND
lines = collect(get_components(Line, c_sys5_pjm_da))
num_line = length(lines)
transformer = collect(get_components(TapTransformer, c_sys5_pjm_da))
num_transformer = length(transformer)
line_loss_ed_ini = zeros(num_line + num_transformer, 1)
for k in 1:num_line
    line_loss_ed_ini[k, 1] = FromTo_Line_ed_ini[k, 3] + ToFrom_Line_ed_ini[k, 3]
end
for k in 1:num_transformer
    line_loss_ed_ini[k+num_line, 1] = FromTo_TapTransformer_ed_ini[k, 3] + ToFrom_TapTransformer_ed_ini[k, 3]
end
E_ed_ini = zeros(num_bus, 1)
for j in 1:num_line
    arc = get_arc(get_component(Line, c_sys5_pjm_da, get_name(lines[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, DF_ed_ini[:, 2]))
    index = findfirst(x -> x == get_name(lines[j]), ToFrom_Line_ed_ini[:, 2])
    E_ed_ini[index1, 1] = E_ed_ini[index1, 1] + line_loss_ed_ini[index, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, DF_ed_ini[:, 2]))
    E_ed_ini[index2, 1] = E_ed_ini[index2, 1] + line_loss_ed_ini[index, 1] / 2
end
for j in 1:num_transformer
    arc = get_arc(get_component(TapTransformer, c_sys5_pjm_da, get_name(transformer[j])))
    bus1 = get_from(arc)
    num1 = get_number(bus1)
    index1 = findfirst(x -> x == num1, parse.(Int, DF_ed_ini[:, 2]))
    index = findfirst(x -> x == get_name(transformer[j]), ToFrom_TapTransformer_ed_ini[:, 2])
    E_ed_ini[index1, 1] = E_ed_ini[index1, 1] + line_loss_ed_ini[index+num_line, 1] / 2
    bus2 = get_to(arc)
    num2 = get_number(bus2)
    index2 = findfirst(x -> x == num2, parse.(Int, DF_ed_ini[:, 2]))
    E_ed_ini[index2, 1] = E_ed_ini[index2, 1] + line_loss_ed_ini[index+num_line, 1] / 2
end
# power injection
injection_ed_ini = value.(ED_models_ini.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])


## run ED with linear loss and fixed HVDC power flow
power_flow_model = ACPowerFlow(;
    exporter=PSSEExportPowerFlow(;
        psse_version=:v33,
        export_dir=psse_export_dir,
        write_comments=true,
        overwrite=true,
    ),
)

template_ed = ProblemTemplate(NetworkModel(PTDFPowerModel; use_slacks=true, duals=[CopperPlateBalanceConstraint], power_flow_evaluation=power_flow_model))
set_device_model!(template_ed, Line, StaticBranchBounds)
set_device_model!(template_ed, TapTransformer, StaticBranchBounds)
set_device_model!(template_ed, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_ed, PowerLoad, StaticPowerLoad)
set_device_model!(template_ed, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_ed, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(
    template_ed,
    DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalDispatch),  # HVDCTwoTerminalLossless
)
set_available!(get_component(HydroDispatch, c_sys5_pjm_da, "201_HYDRO_4"), false)
N_ren=collect(get_components(RenewableNonDispatch, c_sys5_pjm_da))
n_N_ren = length(N_ren)
for j in 1:n_N_ren
    set_available!(N_ren[j], false)
end

ED_models = DecisionModel(
    template_ed,
    c_sys5_pjm_da,
    optimizer=milp_optimizer,
    optimizer_solve_log_print=true,
    name="ED",
    store_variable_names=true,
)

build!(ED_models, output_dir=folder_path)

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
    line_name = idx[1]
    if endswith(line_name, "-double_circuit")
        line_name = replace(line_name, r"-double_circuit$" => "-1")
    end
    line = get_component(Line, c_sys5_pjm_da, line_name)
    from = get_number(get_from(get_arc(line)))
    to = get_number(get_to(get_arc(line)))
    indexx = findfirst(
        x -> x == (from, to),
        axes(ptdf, 2)
    )
    set_normalized_rhs(
        con_tran_ed[idx],
        rhs + sum(
            ptdf[indexx, i] * E_ed_ini[i, 1] / base_power
            for i in 1:num_bus
        )
    )
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
FromTo_Line_ed = aux_variables_ed["PowerFlowBranchActivePowerFromTo__Line"]
ToFrom_Line_ed = aux_variables_ed["PowerFlowBranchActivePowerToFrom__Line"]
FromTo_TapTransformer_ed = aux_variables_ed["PowerFlowBranchActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer_ed = aux_variables_ed["PowerFlowBranchActivePowerToFrom__TapTransformer"]
Ploss_ed = (sum(Matrix(FromTo_Line_ed)[:, 3]) + sum(Matrix(ToFrom_Line_ed)[:, 3]) +
    sum(Matrix(FromTo_TapTransformer_ed)[:, 3]) + sum(Matrix(ToFrom_TapTransformer_ed)[:, 3]))
flow_line_1 = read_variable(results_ed, "FlowActivePowerVariable__Line")
flow_transformer_1 = read_variable(results_ed, "FlowActivePowerVariable__TapTransformer")


# identify cycles
data = PSI.get_power_flow_data(
    PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(ED_models))[1],
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