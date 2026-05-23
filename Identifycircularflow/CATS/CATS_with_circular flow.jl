## There are five parts in the code:
# 1. DCOPF withoutloss
# 2. DCOPF with fixed loss
# 3. DCOPF with quadratic loss
# 4. Identify circular flow based on AC power flow results
# 5. Identify circular flow based on DC power flow results
# Note: The circular flow are identified based on results obtained from 3
# There is circular flow issue

using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using Ipopt
#using Xpress
using Gurobi
using Dates
using JuMP
using PowerFlows
using CSV
using JLD2
using Logging
using InfrastructureSystems
import PowerNetworkMatrices
using XLSX
using JLD2
using Graphs
using SimpleWeightedGraphs
using GraphRecipes
using Plots
using TimeSeries
using Revise
using DataFrames


const PSY = PowerSystems
const PSI = PowerSimulations
const PNM = PowerNetworkMatrices


# Load the saved system #
system = System("CATS_saved_reduced_sys.json")
transform_single_time_series!(
    system,
    Hour(1),  # horizon
    Hour(1),   # interval
)

renewables = collect(get_components(RenewableDispatch, system))
Num_ren = length(renewables)

for row in 1:Num_ren
    new_re_cost = RenewableGenerationCost(;
        variable=CostCurve(LinearCurve(row)),
    )
    set_operation_cost!(get_component(RenewableDispatch, system, get_name(renewables[row])), new_re_cost)
end

hydros = collect(get_components(HydroDispatch, system))
Num_hydro = length(hydros)
for j in 1:Num_hydro
    curve = get_value_curve(get_variable(get_operation_cost(get_component(HydroDispatch, system, get_name(hydros[j])))))
    slope = get_proportional_term(get_function_data(curve)) * 0.1
    constant_term = get_constant_term(get_function_data(curve)) * 0.1

    new_hy_cost = HydroGenerationCost(;
        variable=CostCurve(LinearCurve(-slope, -constant_term)),
        fixed=0,
    )
    set_operation_cost!(get_component(HydroDispatch, system, get_name(hydros[j])), new_hy_cost)
end

thermals = collect(get_components(ThermalStandard, system))
Num_thermals = length(thermals)
for j in 1:Num_thermals
    curve = get_value_curve(get_variable(get_operation_cost(get_component(ThermalStandard, system, get_name(thermals[j])))))
    slope = get_proportional_term(get_function_data(curve)) * 0.01
    constant_term = get_constant_term(get_function_data(curve)) * 0.01

    new_th_cost = ThermalGenerationCost(;
        variable=CostCurve(LinearCurve(-slope, -constant_term)),
        fixed=get_fixed(get_operation_cost(get_component(ThermalStandard, system, get_name(thermals[j])))),
        start_up=get_start_up(get_operation_cost(get_component(ThermalStandard, system, get_name(thermals[j])))),
        shut_down=get_shut_down(get_operation_cost(get_component(ThermalStandard, system, get_name(thermals[j])))),
    )
    set_operation_cost!(get_component(ThermalStandard, system, get_name(thermals[j])), new_th_cost)
end


ptdf_matrix = PTDF(system)

UC_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicUnitCommitment,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

ED_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicUnitCommitment,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    #SynchronousCondenser => SynchronousCondenserBasicDispatch,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)


include("SiennaScripts/build_models.jl")
include("SiennaScripts/utils.jl")
include("SiennaScripts/run_models.jl")
include("SiennaScripts/build_simulations.jl")
include("SiennaScripts/run_simulations.jl")
include("SiennaScripts/add_hvdc.jl")

add_internal_hvdc!(system)
ipopt_optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)
gurobi_global_optimizer = JuMP.optimizer_with_attributes(
    Gurobi.Optimizer,
    "NonConvex" => 2,
    "OutputFlag" => 1,
    "TimeLimit" => 3600,
)


# Use provided PTDF matrices
uc_ptdf_matrix = ptdf_matrix
ed_ptdf_matrix = ptdf_matrix


# Create UC model: determines generator commitment schedules
base_uc_model = make_ptdf_model_without_losses(
    system;
    device_models=UC_MODELS,
    optimizer=gurobi_global_optimizer,
    ptdf=uc_ptdf_matrix,
    name="UC",
    ignore_pf=false,
)

build!(base_uc_model, output_dir=mktempdir())

solve!(base_uc_model)

optimize!(base_uc_model.internal.container.JuMPmodel)
base_uc_objective_value = objective_value(base_uc_model.internal.container.JuMPmodel)

base_uc_bus_injection_values = deepcopy(
    JuMP.value.(
        base_uc_model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
            ActivePowerBalance,
            ACBus,
        }(
            "",
        )]
    ).data,
)
base_uc_results = OptimizationProblemResults(base_uc_model)

base_uc_aux_variables = read_aux_variables(base_uc_results)
base_uc_voltage_values = base_uc_aux_variables["PowerFlowVoltageMagnitude__ACBus"][:, [:value]]
base_uc_voltage_magnitude_df = DataFrame(permutedims(Matrix(base_uc_voltage_values)), Symbol.(1:nrow(base_uc_voltage_values)))

# get loss based on AC power flow calculation
base_uc_line_loss_table = base_uc_aux_variables["PowerFlowBranchActivePowerLoss__Line"]
base_uc_tap_transformer_loss_table = base_uc_aux_variables["PowerFlowBranchActivePowerLoss__Transformer2W"]
base_uc_total_ac_loss = sum(base_uc_line_loss_table[:, 3]) + sum(base_uc_tap_transformer_loss_table[:, 3])
base_uc_injection_values = value.(base_uc_model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
base_uc_total_injection = sum(base_uc_injection_values)

base_uc_thermal_dispatch = read_variable(base_uc_results, "ActivePowerVariable__ThermalStandard")
base_uc_renewable_dispatch = read_variable(base_uc_results, "ActivePowerVariable__RenewableDispatch")
base_uc_hydro_dispatch = read_variable(base_uc_results, "ActivePowerVariable__HydroDispatch")
base_uc_hvdc_flow = read_variable(base_uc_results, "FlowActivePowerVariable__TwoTerminalGenericHVDCLine")



########## run DCOPF with fixed loss ########## 
fixed_loss_uc_model = make_ptdf_model_without_losses(
    system;
    device_models=UC_MODELS,
    optimizer=gurobi_global_optimizer,
    ptdf=uc_ptdf_matrix,
    name="UC_V2",
    ignore_pf=false,
)

build!(fixed_loss_uc_model, output_dir=mktempdir())


# Add loss to copper plate balance
@variable(fixed_loss_uc_model.internal.container.JuMPmodel, fixed_uc_loss_variable <= 0)
ac_buses = collect(get_components(ACBus, system))
slack_bus_positions = PNM.find_slack_positions(ac_buses)
slack_bus_number = get_number(ac_buses[first(slack_bus_positions)])
fixed_uc_balance_constraint = fixed_loss_uc_model.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][slack_bus_number, :]
set_normalized_coefficient(fixed_uc_balance_constraint[1], fixed_uc_loss_variable, 1)

@constraint(fixed_loss_uc_model.internal.container.JuMPmodel,
    fixed_uc_loss_variable == -base_uc_total_ac_loss / 100
)

# Update transmission constraints with fictitious nodal demands
# This ensures branch flow limits account for the additional loading from losses
# by distributing losses to buses and propagating via PTDF
update_transmission_constraints_with_losses!(fixed_loss_uc_model, base_uc_results, system, uc_ptdf_matrix)


solve!(fixed_loss_uc_model)
fixed_loss_uc_results = OptimizationProblemResults(fixed_loss_uc_model)

optimize!(fixed_loss_uc_model.internal.container.JuMPmodel)

fixed_loss_uc_objective_value = objective_value(fixed_loss_uc_model.internal.container.JuMPmodel)

fixed_loss_uc_bus_injection_values = deepcopy(
    JuMP.value.(
        fixed_loss_uc_model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
            ActivePowerBalance,
            ACBus,
        }(
            "",
        )]
    ).data,
)

fixed_uc_aux_variables = read_aux_variables(fixed_loss_uc_results)
fixed_uc_voltage_values = fixed_uc_aux_variables["PowerFlowVoltageMagnitude__ACBus"][:, [:value]]
fixed_uc_voltage_magnitude_df = DataFrame(permutedims(Matrix(fixed_uc_voltage_values)), Symbol.(1:nrow(fixed_uc_voltage_values)))

# get loss based on AC power flow calculation
fixed_uc_line_loss_table = fixed_uc_aux_variables["PowerFlowBranchActivePowerLoss__Line"]
fixed_uc_tap_transformer_loss_table = fixed_uc_aux_variables["PowerFlowBranchActivePowerLoss__Transformer2W"]
fixed_uc_total_ac_loss = sum(fixed_uc_line_loss_table[:, 3]) + sum(fixed_uc_tap_transformer_loss_table[:, 3])
fixed_loss_uc_bus_injection_expr_values = value.(fixed_loss_uc_model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
fixed_uc_total_injection = sum(fixed_loss_uc_bus_injection_values)

fixed_uc_thermal_dispatch = read_variable(fixed_loss_uc_results, "ActivePowerVariable__ThermalStandard")
fixed_uc_renewable_dispatch = read_variable(fixed_loss_uc_results, "ActivePowerVariable__RenewableDispatch")
fixed_uc_hydro_dispatch = read_variable(fixed_loss_uc_results, "ActivePowerVariable__HydroDispatch")
fixed_uc_hvdc_flow = read_variable(fixed_loss_uc_results, "FlowActivePowerVariable__TwoTerminalGenericHVDCLine")



########## run DCOPF with quadratic loss ########## 
# Create ED model: optimizes dispatch given fixed commitments
quadratic_loss_ed_model = make_ptdf_model_without_losses(
    system;
    device_models=ED_MODELS,
    optimizer=ipopt_optimizer,
    ptdf=ed_ptdf_matrix,
    name="ED",
    ignore_pf=false,
)

build!(quadratic_loss_ed_model, output_dir=mktempdir())


# fix binary variables based on the results from fixed_loss_uc_model
fixed_uc_variables = all_variables(fixed_loss_uc_model.internal.container.JuMPmodel)
fixed_uc_solution_by_name = Dict(zip(name.(fixed_uc_variables), value.(fixed_uc_variables)))
ed_variables = all_variables(quadratic_loss_ed_model.internal.container.JuMPmodel)
for ed_var in ed_variables
    if is_binary(ed_var)
        unset_binary(ed_var)
        set_upper_bound(ed_var, fixed_uc_solution_by_name[name(ed_var)])
        set_lower_bound(ed_var, fixed_uc_solution_by_name[name(ed_var)])
    end
end


# Add loss to copper plate balance equation
@variable(quadratic_loss_ed_model.internal.container.JuMPmodel, ed_loss_variable <= 0)
ac_buses = collect(get_components(ACBus, system))
slack_bus_positions = PNM.find_slack_positions(ac_buses)
slack_bus_number = get_number(ac_buses[first(slack_bus_positions)])
ed_balance_constraint = quadratic_loss_ed_model.internal.container.constraints[PowerSimulations.ConstraintKey{CopperPlateBalanceConstraint,System}("")][slack_bus_number, :]
set_normalized_coefficient(ed_balance_constraint[1], ed_loss_variable, 1)

# Extract branch resistance values from system data
branch_resistance, _ = get_RX_vector(system, ptdf_matrix)

# Get current iteration's injection expression (decision variables)
ed_bus_injection_expr = quadratic_loss_ed_model.internal.container.expressions[InfrastructureSystems.Optimization.ExpressionKey{
    ActivePowerBalance,
    ACBus,
}(
    "",
)]
bus_axis = axes(ed_bus_injection_expr, 1)
num_buses = length(bus_axis)


# Extract voltage magnitudes from previous AC power flow solution
bus_voltage_magnitude = get_power_flow_voltage_mag(fixed_loss_uc_results)         # Bus voltages
arc_voltage_magnitude = get_power_flow_arc_voltage_mag(fixed_loss_uc_results, ptdf_matrix)  # Branch endpoint voltages

# select the line for loss calculation
line_loss_name_value = fixed_uc_line_loss_table[:, 2:3]
line_arc_loss_pairs = [
    ((parse(Int, m[1]), parse(Int, m[2])), row.value)
    for row in eachrow(line_loss_name_value)
    for m in [match(r"bus-(\d+)-bus-(\d+)", row.name)]
]
sort!(line_arc_loss_pairs, by=x -> x[2], rev=true)
cumulative_line_losses = cumsum([x[2] for x in line_arc_loss_pairs])
num_selected_lines_80pct = findfirst(x -> x >= 0.9 * cumulative_line_losses[end], cumulative_line_losses)

tap_transformer_loss_name_value = fixed_uc_tap_transformer_loss_table[:, 2:3]
tap_transformer_arc_loss_pairs = [
    ((parse(Int, m[1]), parse(Int, m[2])), row.value)
    for row in eachrow(tap_transformer_loss_name_value)
    for m in [match(r"bus-(\d+)-bus-(\d+)", row.name)]
]
sort!(tap_transformer_arc_loss_pairs, by=x -> x[2], rev=true)
cumulative_tap_transformer_losses = cumsum([x[2] for x in tap_transformer_arc_loss_pairs])
num_selected_tap_transformers_80pct = findfirst(x -> x >= 0.9 * cumulative_tap_transformer_losses[end], cumulative_tap_transformer_losses)

selected_line_arcs = [x[1] for x in line_arc_loss_pairs[1:num_selected_lines_80pct]]
selected_tap_transformer_arcs = [x[1] for x in tap_transformer_arc_loss_pairs[1:num_selected_tap_transformers_80pct]]
selected_loss_arcs = vcat(selected_line_arcs, selected_tap_transformer_arcs)

ptdf_arc_axis = ptdf_matrix.axes[2]
selected_arc_indices = [findfirst(x -> x == pair, ptdf_arc_axis) for pair in selected_loss_arcs]
selected_branch_resistance = branch_resistance[selected_arc_indices]
selected_arc_voltage_magnitude = arc_voltage_magnitude[selected_arc_indices]
ptdf_data_matrix = PowerNetworkMatrices.get_ptdf_data(ptdf_matrix)
selected_ptdf_rows = ptdf_data_matrix[selected_arc_indices, :]
num_selected_arcs = length(selected_arc_indices)


# Create quadratic loss constraint for each time period
@variable(quadratic_loss_ed_model.internal.container.JuMPmodel, selected_arc_flow_for_loss[1:num_selected_arcs])
@constraint(quadratic_loss_ed_model.internal.container.JuMPmodel, [k in 1:num_selected_arcs],
    selected_arc_flow_for_loss[k] ==
    sum(
        (selected_arc_voltage_magnitude[k] / bus_voltage_magnitude[j]) *
        selected_ptdf_rows[k, j] *
        ed_bus_injection_expr[bus_axis[j], 1]
        for j in 1:num_buses
        if abs(selected_ptdf_rows[k, j]) > 1e-8
    )
)

@constraint(quadratic_loss_ed_model.internal.container.JuMPmodel,
    ed_loss_variable + (sum(selected_branch_resistance[k] * selected_arc_flow_for_loss[k]^2 for k in 1:num_selected_arcs)) == 0
)

# Update transmission constraints with fictitious nodal demands
# This ensures branch flow limits account for the additional loading from losses
# by distributing losses to buses and propagating via PTDF
update_transmission_constraints_with_losses!(quadratic_loss_ed_model, base_uc_results, system, ed_ptdf_matrix)

solve!(quadratic_loss_ed_model)
quadratic_loss_ed_results = OptimizationProblemResults(quadratic_loss_ed_model)

optimize!(quadratic_loss_ed_model.internal.container.JuMPmodel)

quadratic_loss_ed_objective_value = objective_value(quadratic_loss_ed_model.internal.container.JuMPmodel)

# obtain the injection at each bus
ed_thermal_dispatch = read_variable(quadratic_loss_ed_results, "ActivePowerVariable__ThermalStandard")
ed_renewable_dispatch = read_variable(quadratic_loss_ed_results, "ActivePowerVariable__RenewableDispatch")
ed_hydro_dispatch = read_variable(quadratic_loss_ed_results, "ActivePowerVariable__HydroDispatch")
ed_hvdc_flow = read_variable(quadratic_loss_ed_results, "FlowActivePowerVariable__TwoTerminalGenericHVDCLine")
renewable_dispatch_limit = read_parameter(quadratic_loss_ed_results, "ActivePowerTimeSeriesParameter__RenewableDispatch")
hydro_dispatch_limit = read_parameter(quadratic_loss_ed_results, "ActivePowerTimeSeriesParameter__HydroDispatch")
load_active_power_param = read_parameter(quadratic_loss_ed_results, "ActivePowerTimeSeriesParameter__PowerLoad")
read_variable(quadratic_loss_ed_results, "SystemBalanceSlackDown__System")

ed_aux_variables = read_aux_variables(quadratic_loss_ed_results)
ed_line_ac_loss_vector = ed_aux_variables["PowerFlowBranchActivePowerLoss__Line"][:, 3]
ed_tap_transformer_ac_loss_vector = ed_aux_variables["PowerFlowBranchActivePowerLoss__Transformer2W"][:, 3]
ed_total_ac_loss = sum(ed_line_ac_loss_vector) + sum(ed_tap_transformer_ac_loss_vector)
ed_bus_injection_values = value.(quadratic_loss_ed_model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{ActivePowerBalance,ACBus}("")])
ed_total_modeled_loss_mw = sum(ed_bus_injection_values) * 100

## identify circular flow based on AC power flow
include("SiennaScripts/mapped_indices.jl")
include("SiennaScripts/circular_flows.jl")
data = PSI.get_power_flow_data(
    only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(quadratic_loss_ed_model))),
)
base_power = get_base_power(system)

#res_vars = read_variables(res_old_uc)
time_step1 = 1
G = build_graph(data; time_step=time_step1)
time_step2 = DateTime("2019-01-01T12:00:00")
add_hvdc_edges!(G, system, ed_hvdc_flow, data; time_step=time_step2)
branches = collect(get_components(ACBranch, system))
C = find_circular_flows(G, data, branches)
num_cycle = length(C)
circular_flow = zeros(num_cycle, 1)
for j in 1:num_cycle
    circular_flow[j, 1] = minimum(C[j, 1].branch_flows)
end

XLSX.openxlsx("circular flow.xlsx", mode="w") do xf
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "cycles")
    XLSX.writetable!(sheet1, Tables.table(Matrix(circular_flow)))
end


## identify circular flow based on DC power flow results
ed_bus_injection_values = deepcopy(
    JuMP.value.(
        quadratic_loss_ed_model.internal.container.expressions[PSY.InfrastructureSystems.Optimization.ExpressionKey{
            ActivePowerBalance,
            ACBus,
        }(
            "",
        )]
    ).data,
)
# need to recalculate the power flow, because only several line's flow is saved in DCOPF
num_lines = length(ptdf_matrix.axes[2])
line_flow_ed = zeros(num_lines, 1)
for j in 1:num_lines
    line_flow_ed[j, 1] = sum(ptdf_matrix[j, k] * ed_bus_injection_values[k] for k in 1:num_buses)
end
from = Vector{Int64}()
to = Vector{Int64}()
w = Vector{Float64}()

line_names = ptdf_matrix.axes[2]
lines = collect(get_components(Line, system))

line_name_set = [
    ((parse(Int, m[1]), parse(Int, m[2])))
    for row in eachrow(line_loss_name_value)
    for m in [match(r"bus-(\d+)-bus-(\d+)", row.name)]
]
transformers_name_set = [
    ((parse(Int, m[1]), parse(Int, m[2])))
    for row in eachrow(tap_transformer_loss_name_value)
    for m in [match(r"bus-(\d+)-bus-(\d+)", row.name)]
]

num_line = length(line_names)
for j in 1:num_line
    line_name = line_names[j]
    index1 = findfirst(x -> x == line_name, line_name_set)
    if !isnothing(index1)
        arc = get_arc(get_component(Line, system, line_loss_name_value[index1, 1]))
        bus1 = get_number(get_from(arc))
        bus2 = get_number(get_to(arc))
        if line_flow_ed[j, 1] > 0
            push!(from, bus1)
            push!(to, bus2)
            push!(w, line_flow_ed[j, 1])
        else
            push!(from, bus2)
            push!(to, bus1)
            push!(w, -line_flow_ed[j, 1])
        end
    else
        index2 = findfirst(x -> x == line_name, transformers_name_set)
        arc = get_arc(get_component(Transformer2W, system, tap_transformer_loss_name_value[index2, 1]))
        bus1 = get_number(get_from(arc))
        bus2 = get_number(get_to(arc))
        if line_flow_ed[j, 1] > 0
            push!(from, bus1)
            push!(to, bus2)
            push!(w, line_flow_ed[j, 1])
        else
            push!(from, bus2)
            push!(to, bus1)
            push!(w, -line_flow_ed[j, 1])
        end
    end
end


if ed_hvdc_flow[1, 3] > 0
    push!(from, 8335)
    push!(to, 8814)
    push!(w, ed_hvdc_flow[1, 3])
else
    push!(from, 8814)
    push!(to, 8335)
    push!(w, -ed_hvdc_flow[1, 3])
end

if ed_hvdc_flow[2, 3] > 0
    push!(from, 1819)
    push!(to, 1258)
    push!(w, ed_hvdc_flow[2, 3])
else
    push!(from, 1258)
    push!(to, 1819)
    push!(w, -ed_hvdc_flow[2, 3])
end

G = SimpleWeightedDiGraph(from, to, w)

cycles = simplecycles(G)

