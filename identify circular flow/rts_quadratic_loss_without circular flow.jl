# solve the SCUC problem considering quadratic loss for RTS
# consider renewable generators
# consider HVDC
# code to get the results of cost-minimization model
# no congestion
# there is no circular flow issue


using Pkg
Pkg.activate(@__DIR__)
using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using PowerFlows
using TimeSeries
#using Ipopt
using Gurobi
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
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(1), Hour(1))

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

set_reactive_power_limits_from!(hvdc1, (min=-100, max=100))
set_reactive_power_limits_to!(hvdc1, (min=-100, max=100))
set_active_power_limits_from!(hvdc1, (min=-100, max=100))
set_active_power_limits_to!(hvdc1, (min=-100, max=100))

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

milp_optimizer = optimizer_with_attributes(
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

models = DecisionModel(
    template_uc,
    c_sys5_pjm_da,
    optimizer=milp_optimizer,
    optimizer_solve_log_print=true,
    name="UC",
    store_variable_names=true,
)

build!(models, output_dir=mktempdir())


@variable(models.internal.container.JuMPmodel, loss <= 0)
buses = collect(get_components(ACBus, c_sys5_pjm_da))
indexx = PNM.find_slack_positions(buses)
ref_bus = get_number(buses[first(indexx)])
con_bal = models.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][ref_bus, :]
set_normalized_coefficient(con_bal[1], loss, 1)
injection = models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")]
AX = axes(injection, 1)
ptdf = PTDF(c_sys5_pjm_da)
AA = axes(ptdf, 2)
lines = collect(get_components(Line, c_sys5_pjm_da))
num_line = length(lines)
transformer = collect(get_components(TapTransformer, c_sys5_pjm_da))
num_transformer = length(transformer)
R = zeros(num_line + num_transformer, 1)
for j in 1:num_line
    index1 = findfirst(x -> x == get_name(lines[j]), AA)
    if index1 !== nothing
        R[index1, 1] = get_r(lines[j])
    end
end
for j in 1:num_transformer
    index1 = findfirst(x -> x == get_name(transformer[j]), AA)
    if index1 !== nothing
        R[index1, 1] = get_r(transformer[j])
    end
end
Num_Jen = length(AX)

@constraint(models.internal.container.JuMPmodel, loss == -sum(R[k] * (sum(ptdf[k, j] * injection[AX[j], 1] for j in 1:Num_Jen))^2 for k in 1:num_line+num_transformer))

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

loss_0 = sum(injection_1)

aux_variables = read_aux_variables(results_uc_1)
# get loss factor based on AC power flow calculation
DF = aux_variables["PowerFlowLossFactors__ACBus"] # penalty factor
buses = collect(get_components(ACBus, c_sys5_pjm_da))
num_bus = length(buses)
LF = ones(num_bus, 1) + Matrix(DF)[1, 2:end]
indexx = PNM.find_slack_positions(buses)
ref_bus = get_number(buses[first(indexx)])
indexx = findfirst(x -> x == ref_bus, parse.(Int, names(DF)[2:end]))
LF[indexx] = 0

# get estimated power losses based on AC power flow calculation
FromTo_Line = aux_variables["PowerFlowLineActivePowerFromTo__Line"]
ToFrom_Line = aux_variables["PowerFlowLineActivePowerToFrom__Line"]
FromTo_TapTransformer = aux_variables["PowerFlowLineActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer = aux_variables["PowerFlowLineActivePowerToFrom__TapTransformer"]
Ploss_est = (sum(Matrix(FromTo_Line)[1, 2:end]) + sum(Matrix(ToFrom_Line)[1, 2:end]) +
             sum(Matrix(FromTo_TapTransformer)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer)[1, 2:end]))


# get line flow from DCOPF calculation
flow_line_1 = read_variable(results_uc_1, "FlowActivePowerVariable__Line")
flow_transformer_1 = read_variable(results_uc_1, "FlowActivePowerVariable__TapTransformer")


HVDC_flow_1 = read_variable(results_uc_1, "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine")
HVDC_loss = read_variable(results_uc_1, "HVDCLosses__TwoTerminalGenericHVDCLine")



## identify circular flow based on AC power flow results
data = PSI.get_power_flow_data(
    only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(models))),
)
base_power = get_base_power(c_sys5_pjm_da)

res_vars = read_variables(results_uc_1)
time_step = 1
G = build_graph(data; time_step=time_step)
add_hvdc_edges!(G, c_sys5_pjm_da, res_vars, data; time_step=time_step)
branches = collect(get_components(ACBranch, c_sys5_pjm_da))
C = find_circular_flows(G, data, branches)
num_cycle = length(C)
circular_flow = zeros(num_cycle, 1)
for j in 1:num_cycle
    circular_flow[j, 1] = minimum(C[j, 1].branch_flows)
end


## ideitify circular flow based on DCOPF results
flow_line_3 = read_variable(results_uc_1, "FlowActivePowerVariable__Line")
flow_transformer_3 = read_variable(results_uc_1, "FlowActivePowerVariable__TapTransformer")
HVDC_flow_3 = read_variable(results_uc_1, "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine")
from = Vector{Int64}()
to = Vector{Int64}()
w = Vector{Float64}()

line_names = names(flow_line_3)
for j in 1:num_line
    line_name = get_name(lines[j])
    index1 = findfirst(x -> x == line_name, line_names)
    arc = get_arc(get_component(Line, c_sys5_pjm_da, get_name(lines[j])))
    bus1 = get_number(get_from(arc))
    bus2 = get_number(get_to(arc))
    if flow_line_3[1, index1] > 0
        push!(from, bus1)
        push!(to, bus2)
        push!(w, flow_line_3[1, index1])
    else
        push!(from, bus2)
        push!(to, bus1)
        push!(w, -flow_line_3[1, index1])
    end
end

tran_names = names(flow_transformer_3)
for j in 1:num_transformer
    tran_name = get_name(transformer[j])
    index2 = findfirst(x -> x == tran_name, tran_names)
    arc = get_arc(get_component(TapTransformer, c_sys5_pjm_da, get_name(transformer[j])))
    bus1 = get_number(get_from(arc))
    bus2 = get_number(get_to(arc))
    if flow_transformer_3[1, index2] > 0
        push!(from, bus1)
        push!(to, bus2)
        push!(w, flow_transformer_3[1, index2])
    else
        push!(from, bus2)
        push!(to, bus1)
        push!(w, -flow_transformer_3[1, index2])
    end
end

if HVDC_flow_3[1, 2] > 0
    push!(from, 113)
    push!(to, 316)
    push!(w, HVDC_flow_3[1, 2])
else
    push!(from, 316)
    push!(to, 113)
    push!(w, -HVDC_flow_3[1, 2])
end

G = SimpleWeightedDiGraph(from, to, w)

cycles = simplecycles(G)