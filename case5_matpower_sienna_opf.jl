using Pkg
Pkg.activate(".")
using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using Ipopt
using PowerFlows
using TimeSeries

const PSY = PowerSystems

sys = build_system(PSISystems, "c_sys5_pjm")

# Update Data to Match Matpower case
for x in collect(get_components(RenewableDispatch, sys))
    remove_component!(sys, x)
end
for x in get_components(ACBus, sys)
    set_voltage_limits!(x, (min = 0.9, max = 1.1))
end
for x in get_components(PowerLoad, sys)
    set_max_active_power!(x, get_active_power(x))
end
set_active_power!(get_component(ThermalStandard, sys, "Brighton"), 4.6651)
set_active_power!(get_component(ThermalStandard, sys, "Solitude"), 3.2349)

# Update Load if wanted
l = first(get_components(PowerLoad, sys))
ts = get_time_series_array(SingleTimeSeries, l, "max_active_power")
timestamps_vals = timestamp(ts)
load_values = ones(168)
load_timearray = TimeArray(timestamps_vals, load_values)

new_load_timeseries = SingleTimeSeries(;
    name = "max_active_power",
    data = load_timearray,
    scaling_factor_multiplier = get_max_active_power)

#remove_time_series!(sys, SingleTimeSeries)
#loads = collect(get_components(PowerLoad, sys))
#add_time_series!(sys, loads, new_load_timeseries)


PSY.transform_single_time_series!(sys, Hour(24), Hour(24))

# Set ACP Template
template_ed = ProblemTemplate(NetworkModel(ACPPowerModel))
set_device_model!(template_ed, ThermalStandard, ThermalDispatchNoMin)
set_device_model!(template_ed, Line, StaticBranchUnbounded)
set_device_model!(template_ed, PowerLoad, StaticPowerLoad)


nlp_optimizer = optimizer_with_attributes(Ipopt.Optimizer)
model = DecisionModel(
    template_ed,
    sys,
    optimizer=nlp_optimizer,
    name="ED",
    system_to_file=false,
    initialize_model=false,
)
build!(model, output_dir = mktempdir())

solve!(model)

res = OptimizationProblemResults(model)

vmag = read_variable(res, "VoltageMagnitude__ACBus")
th_power = read_variable(res, "ActivePowerVariable__ThermalStandard")
th_rpower = read_variable(res, "ReactivePowerVariable__ThermalStandard")
load_param = read_parameter(res, "ActivePowerTimeSeriesParameter__PowerLoad")
load_param_re = read_parameter(res, "ReactivePowerTimeSeriesParameter__PowerLoad")

show_components(sys, ThermalStandard, [:bus])