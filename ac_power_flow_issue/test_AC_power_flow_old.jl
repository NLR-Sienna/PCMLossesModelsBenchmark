##### Cost-minimization problem
# UC formulated as DCOPF with fixed HVDC power flows (HVDC_value) and fixed losses (Ploss_ori);
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
# power injection
injection_uc = value.(UC_models.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
DC_loss = sum(injection_uc[:, 1])*base_power

aux_variables_uc = read_aux_variables(results_uc)
# estimated losses
FromTo_Line = aux_variables_uc["PowerFlowLineActivePowerFromTo__Line"]
ToFrom_Line = aux_variables_uc["PowerFlowLineActivePowerToFrom__Line"]
FromTo_TapTransformer = aux_variables_uc["PowerFlowLineActivePowerFromTo__TapTransformer"]
ToFrom_TapTransformer = aux_variables_uc["PowerFlowLineActivePowerToFrom__TapTransformer"]
AC_loss = (sum(Matrix(FromTo_Line)[1, 2:end]) + sum(Matrix(ToFrom_Line)[1, 2:end]) +
    sum(Matrix(FromTo_TapTransformer)[1, 2:end]) + sum(Matrix(ToFrom_TapTransformer)[1, 2:end]))
# loss factor
DF = aux_variables_uc["PowerFlowLossFactors__ACBus"] # penalty factor
LF = ones(num_bus, 1) + Matrix(DF)[1, 2:end]
indexx_ref = findfirst(x -> x == ref_bus, parse.(Int, names(DF)[2:end]))
LF[indexx_ref] = 0

println("DC_loss: $DC_loss")
println("AC_loss: $AC_loss")
