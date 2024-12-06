using Pkg
Pkg.activate(".")
using Revise
using PowerSimulations
using PowerSystems
using PowerSystemCaseBuilder
const PSY = PowerSystems
const PSI = PowerSimulations
const PSB = PowerSystemCaseBuilder
using HiGHS
using Ipopt
using Logging
using Dates
using CSV
using DataFrames
#using PowerGraphics
# plotlyjs()

c_sys5_pjm_da = PSB.build_system(PSISystems, "c_sys5_pjm")
PSY.transform_single_time_series!(c_sys5_pjm_da, Hour(24), Hour(24))
c_sys5_pjm_rt = PSB.build_system(PSISystems, "c_sys5_pjm_rt")
PSY.transform_single_time_series!(c_sys5_pjm_rt, Hour(2), Hour(1))

for sys in [c_sys5_pjm_da, c_sys5_pjm_rt]
    th = get_component(ThermalStandard, sys, "Park City")
    set_active_power_limits!(th, (min = 0.1, max = 1.7))
    set_status!(th, false)
    set_active_power!(th, 0.0)
    c = get_operation_cost(th)
    set_start_up!(c, 1500.0)
    set_shut_down!(c, 75.0)
    set_time_at_status!(th, 1)

    th = get_component(ThermalStandard, sys, "Alta")
    set_time_limits!(th, (up = 5, down = 1))
    set_active_power_limits!(th, (min = 0.05, max = 0.4))
    set_active_power!(th, 0.05)
    c = get_operation_cost(th)
    set_start_up!(c, 400.0)
    set_shut_down!(c, 200.0)
    set_time_at_status!(th, 2)

    th = get_component(ThermalStandard, sys, "Brighton")
    set_active_power_limits!(th, (min = 2.0, max = 6.0))
    c = get_operation_cost(th)
    set_active_power!(th, 4.88041)
    set_start_up!(c, 5000.0)
    set_shut_down!(c, 3000.0)

    th = get_component(ThermalStandard, sys, "Sundance")
    set_active_power_limits!(th, (min = 1.0, max = 2.0))
    set_time_limits!(th, (up = 5, down = 1))
    set_active_power!(th, 2.0)
    c = get_operation_cost(th)
    set_start_up!(c, 4000.0)
    set_shut_down!(c, 2000.0)
    set_time_at_status!(th, 1)

    th = get_component(ThermalStandard, sys, "Solitude")
    set_active_power_limits!(th, (min = 1.0, max = 5.2))
    set_ramp_limits!(th, (up = 0.0052, down = 0.0052))
    set_active_power!(th, 2.0)
    c = get_operation_cost(th)
    set_start_up!(c, 3000.0)
    set_shut_down!(c, 1500.0)
end

milp_optimizer = optimizer_with_attributes(HiGHS.Optimizer)
nlp_optimizer = optimizer_with_attributes(Ipopt.Optimizer)

template_uc = template_unit_commitment()
set_network_model!(
    template_uc,
    NetworkModel(
        PTDFPowerModel;
    ),
)

set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
template_ed = deepcopy(template_uc)
set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_network_model!(
    template_ed,
    NetworkModel(
        ACPPowerModel,
        use_slacks = true
    ),
)

models = SimulationModels(
    decision_models=[
        DecisionModel(
            template_uc,
            c_sys5_pjm_da,
            optimizer=milp_optimizer,
            name="UC",
            system_to_file=true,
        ),
        DecisionModel(
            template_ed,
            c_sys5_pjm_rt,
            optimizer=nlp_optimizer,
            name="ED",
            system_to_file=true,
            initialize_model=false,
        ),
    ],
)

sequence = SimulationSequence(
    models=models,
    feedforwards=Dict(
        "ED" => [
            SemiContinuousFeedforward(
                component_type=ThermalStandard,
                source=OnVariable,
                affected_values=[ActivePowerVariable, ReactivePowerVariable],
            ),
        ],
    ),
    ini_cond_chronology=InterProblemChronology(),
)

sim = Simulation(
    name="5Bus_sim",
    steps=6,
    models=models,
    sequence=sequence,
    simulation_folder=mktempdir(),
)

build!(sim; console_level=Logging.Info)
execute!(sim)

res_sim = SimulationResults(sim)
results_uc = get_decision_problem_results(res_sim, "UC")
results_ed = get_decision_problem_results(res_sim, "ED")

plot_fuel(results_uc; stair = true)
plot_fuel(results_ed)

f_var_uc = read_realized_variable(results_uc, "FlowActivePowerVariable__Line")
f_var_ed = read_realized_variable(results_ed, "FlowActivePowerFromToVariable__Line")

plot_dataframe(f_var_uc; stair = true)
plot_dataframe(f_var_ed)

p_re_var = read_realized_variable(results_uc, ActivePowerVariable, RenewableDispatch)
ren_data_uc = read_realized_parameter(results_uc, "ActivePowerTimeSeriesParameter__RenewableDispatch")
ren_data_ed = read_realized_variable(results_ed, "SystemBalanceSlackUp__Bus__P")
ren_data_ed = read_realized_variable(results_ed, "VoltageMagnitude__Bus")
plot_dataframe(ren_data_ed)


ren_data .- p_re_var



plexos_results_folder = "PLEXOS-5-BUS/Base_csv_results"
p_var_plexos = CSV.read(joinpath(plexos_results_folder, "generation.csv"), DataFrame)
rename!(p_var_plexos, Dict(:timestamp => :DateTime))
rename!(p_var_plexos, Dict(:Park_City => Symbol("Park City")))
rename!(p_var_plexos, Dict(:SolarBusC => :PVBus5))
rename!(p_var_plexos, Dict(:WindBusA => :WindBus1))

thermal_cols = ["DateTime", "Alta", "Brighton", "Park City", "Solitude", "Sundance"]
thermal_gen_plexos = p_var_plexos[!, thermal_cols]

var_cost_exp = read_realized_expression(results_uc, "ProductionCostExpression__ThermalStandard")
plexos_var_cost = deepcopy(thermal_gen_plexos)
plexos_var_cost[!, "Alta"] = plexos_var_cost[!, "Alta"]*14.0
plexos_var_cost[!, "Brighton"] = plexos_var_cost[!, "Brighton"]*10.0
plexos_var_cost[!, "Sundance"] = plexos_var_cost[!, "Sundance"]*40.0
plexos_var_cost[!, "Park City"] = plexos_var_cost[!, "Park City"]*15.0
plexos_var_cost[!, "Solitude"] = plexos_var_cost[!, "Solitude"]*30.0

cost_diff = sum(eachcol(plexos_var_cost)[2:end]) .- sum(eachcol(var_cost_exp)[2:end])
plot(p_th_var[!, 1], cost_diff)

start_th_var = read_realized_variable(results_uc, StartVariable, ThermalStandard)
plot_dataframe(start_th_var)
off_th_var = read_realized_variable(results_uc, StopVariable, ThermalStandard)
plot_dataframe(off_th_var)

ren_cols = ["DateTime", "PVBus5", "WindBus1"]
ren_gen_plexos = p_var_plexos[!, ren_cols]

line_flow_plexos = CSV.read(joinpath(plexos_results_folder, "flow.csv"), DataFrame)
rename!(line_flow_plexos, names(f_var))

difference_p_th = thermal_gen_plexos[!, 2:end] .- p_th_var[!, 2:end]
plot_dataframe(difference_p_th, p_th_var[!, 1])

diff_ren = ren_gen_plexos[!, 2:end] .- p_re_var[!, 2:end]
plot_dataframe(diff_ren, p_th_var[!, 1])

plot_fuel(results_uc; stair = true)

wind_plexos = CSV.read("PLEXOS-5-BUS/5_bus_inputs/da_wind.csv", DataFrame)
solar_plexos = CSV.read("PLEXOS-5-BUS/5_bus_inputs/da_solar.csv", DataFrame)
load_plexos = CSV.read("PLEXOS-5-BUS/5_bus_inputs/da_load.csv", DataFrame)
load_data = read_realized_parameter(results_uc, "ActivePowerTimeSeriesParameter__PowerLoad")

maximum(ren_data[!, "WindBus1"] .- wind_plexos[!, "WindBusA"])
maximum(ren_data[!, "PVBus5"] .- solar_plexos[!, "SolarBusC"])
maximum(-1*sum(eachcol(load_data)[2:end]) .- load_plexos[!, "region1"])

plot_dataframe(p_th_var; stair = true)
plot_dataframe(p_re_var; stair = true)
plot_dataframe(f_var; stair = true)


on_th_var = read_realized_variable(results_uc, OnVariable, ThermalStandard)
on_th_var[1:26, :]

plot_dataframe(on_th_var)
