using Graphs
using SimpleWeightedGraphs
using GraphRecipes
using Plots

using HiGHS
using Dates
using TimeSeries

import Xpress
import Ipopt

import HSL_jll # works only after getting the HSL license and the files for HSL_jll; must install: ] dev <path_to_HSL_jll>
using HSL # works only after getting the HSL license and the files for HSL_jll (]dev HSL_jll first)

import JuMP


using PowerSystems
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using PowerFlows
using PowerSimulations
using HydroPowerSimulations
PSI = PowerSimulations

const PSY = PowerSystems

@assert LIBHSL_isfunctional()

"find cycles in a small system with HVDC"

_check_name(sys, name, component) = name  # placeholder for now

"""    
    _add_simple_bus!(sys::System, number::Int, bus_type::ACBusTypes, base_voltage::Number, voltage_magnitude::Float64=1.0, voltage_angle::Float64=0.0)
    Simplified function to create and add a bus to the system with the given parameters.
"""
function _add_simple_bus!(
    sys::System,
    number::Int,
    bus_type::ACBusTypes,
    base_voltage::Number,
    voltage_magnitude::Float64 = 1.0,
    voltage_angle::Float64 = 0.0,
)
    bus = ACBus(;
        number = number,
        name = _check_name(sys, "bus_$number", ACBus),
        bustype = bus_type,
        angle = voltage_angle,
        magnitude = voltage_magnitude,
        voltage_limits = (0.8, 1.2),
        base_voltage = Float64(base_voltage),
    )
    add_component!(sys, bus)
    return bus
end


"""    
    _add_simple_load!(sys::System, bus::ACBus, active_power::Number, reactive_power::Number)
    Simplified function to create and add a load to the system with the given parameters.
"""
function _add_simple_load!(
    sys::System,
    bus::ACBus,
    active_power::Number,
    reactive_power::Number,
    base_power::Number = 100.0,
)
    load = PowerLoad(;
        name = _check_name(sys, "load_$(get_number(bus))", PowerLoad),
        available = true,
        bus = bus,
        active_power = Float64(active_power), # Per-unitized by device base_power
        reactive_power = Float64(reactive_power), # Per-unitized by device base_power
        base_power = Float64(base_power), # MVA
        max_active_power = 1.0, # 10 MW per-unitized by device base_power
        max_reactive_power = 1.0,
    )

    add_component!(sys, load)
    return load
end


"""    
    _add_simple_thermal_standard!(sys::System, bus::ACBus, active_power::Number=0.0, reactive_power::Number=0.0)
    Simplified function to create and add a thermal standard generator to the system with the given parameters.
"""
function _add_simple_thermal_standard!(
    sys::System,
    bus::ACBus,
    active_power::Number,
    reactive_power::Number,
)
    gen = ThermalStandard(;
        name = _check_name(sys, "thermal_standard_$(get_number(bus))", ThermalStandard),
        available = true,
        status = true,
        bus = bus,
        active_power = Float64(active_power),
        reactive_power = Float64(reactive_power),
        rating = 1.0,
        active_power_limits = (min = 0, max = 100),
        reactive_power_limits = (min = -100, max = 100),
        ramp_limits = nothing,
        operation_cost = ThermalGenerationCost(CostCurve(LinearCurve(10)), 0., 0., 0.),
        base_power = 100.0,
        time_limits = nothing,
        prime_mover_type = PrimeMovers.OT,
        fuel = ThermalFuels.OTHER,
        services = Device[],
        dynamic_injector = nothing,
        ext = Dict{String, Any}(),
    )
    add_component!(sys, gen)
    return gen
end

"""    
    _add_simple_line!(sys::System, bus1::ACBus, bus2::ACBus, r::Float64=1e-3, x::Float64=1e-3, b::Float64=0.0)
    Simplified function to create and add a line to the system with the given parameters.
"""
function _add_simple_line!(
    sys::System,
    bus1::ACBus,
    bus2::ACBus,
    r::Float64 = 1e-3,
    x::Float64 = 1e-3,
    b::Float64 = 0.0,
)
    line = Line(;
        name = _check_name(sys, "line_$(get_number(bus1))_$(get_number(bus2))", Line),
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = bus1, to = bus2),
        r = r,
        x = x,
        b = (from = b / 2, to = b / 2),
        rating = 2.0,
        angle_limits = (min = -pi / 2, max = pi / 2),
    )
    add_component!(sys, line)
    return line
end


function minimal_system()
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PQ, 230, 1.0, 0.0)
    b3 = _add_simple_bus!(sys, 3, ACBusTypes.PQ, 230, 1.0, 0.0)
    b4 = _add_simple_bus!(sys, 4, ACBusTypes.PQ, 230, 1.0, 0.0)

    _add_simple_thermal_standard!(sys, b1, 0, 0)

    l1 = _add_simple_line!(sys, b1, b2)
    l2 = _add_simple_line!(sys, b2, b3)
    l3 = _add_simple_line!(sys, b2, b4)

    hvdc1 = TwoTerminalGenericHVDCLine(
        name = "HVDC1",
        available = true,
        active_power_flow = 0.0,
        arc=Arc(from=b3, to=b4),
        active_power_limits_from=(min=-10.0, max=10.0),
        active_power_limits_to=(min=-10.0, max=10.0),
        reactive_power_limits_from=(min=-10.0, max=10.0),
        reactive_power_limits_to=(min=-10.0, max=10.0),
        loss=LinearCurve(0.),
    )
    
    add_component!(sys, hvdc1)

    g1 = RenewableDispatch(
        name = "Wind1",
        available = true,
        bus=b3,
        active_power=0.0,
        reactive_power=0.0,
        rating=1.0,
        prime_mover_type=PrimeMovers.WT,
        reactive_power_limits=(min=-1., max=1.),
        power_factor=0.95,
        operation_cost=RenewableGenerationCost(
            curtailment_cost=CostCurve(LinearCurve(100.0)),
            variable=CostCurve(LinearCurve(-10.0)),
            fixed=0.0,
        ),
        base_power=100.0,
    )
    add_component!(sys, g1)

    horizon=24

    ts_data = [1.0, 0.99, 0.99, 1.0, 0.99, 0.99, 0.99, 0.98, 0.95, 0.92, 0.90, 0.88, 0.84, 0.76,
           0.65, 0.52, 0.39, 0.28, 0.19, 0.15, 0.13, 0.11, 0.09, 0.06,]
    time_stamps = range(DateTime("2020-01-01"); step = Hour(1), length = horizon)
    time_series_data_raw = TimeArray(time_stamps, ts_data)
    time_series = SingleTimeSeries(; name = "max_active_power", data = time_series_data_raw)

    add_time_series!(sys, g1, time_series)

    ld1 = _add_simple_load!(sys, b4, 0.5, 0.25, 20.)

    time_series = SingleTimeSeries(; name = "max_active_power", data = time_series_data_raw)

    add_time_series!(sys, ld1, time_series)

    transform_single_time_series!(sys, Hour(horizon), Hour(horizon))

    return sys
end

sys_uc = minimal_system()
sys_ed = minimal_system()

PSY.transform_single_time_series!(sys_uc, Hour(2), Hour(2))
PSY.transform_single_time_series!(sys_ed, Hour(1), Hour(1))

template_uc = ProblemTemplate()

set_device_model!(template_uc, Line, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)
set_device_model!(template_uc, TapTransformer, StaticBranch)

set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)  # RenewableFullDispatch
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_uc, RenewableNonDispatch, FixedOutput)

set_device_model!(
        template_uc,
        DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalDispatch),  # HVDCTwoTerminalLossless
    )

# set_service_model!(template_uc, VariableReserve{ReserveUp}, RangeReserve)
# set_service_model!(template_uc, VariableReserve{ReserveDown}, RangeReserve)

set_network_model!(template_uc, NetworkModel(DCPPowerModel; use_slacks=false))

solver_highs = optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.5)

solver_xpress = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.02)

solver_ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    "print_level" => 1,
    "hsllib" => HSL_jll.libhsl_path, # uncomment after getting the HSL license and the files for HSL_jll
    "linear_solver" => "ma57", # uncomment after getting the HSL license and the files for HSL_jll
    "tol" => 1e-6,
    "acceptable_tol" => 1e-5,
    )

problem_uc = DecisionModel(template_uc, sys_uc; optimizer = solver_xpress, name="UC")


##########
network_model_ed = NetworkModel(
        ACPPowerModel; 
        # DCPPowerModel;
        use_slacks=false, 
        power_flow_evaluation=PowerFlows.ACPowerFlow(
            ;
            # calculate_loss_factors=true, 
            # check_reactive_power_limits=q_lim,
            # generator_slack_participation_factors=ds ? Dict(get_name(x) => 1.0 for x in get_components(Generator, system)) : nothing,
        ),
    )

template_ed = ProblemTemplate(network_model_ed)

set_device_model!(template_ed, Line, StaticBranch)
set_device_model!(template_ed, Transformer2W, StaticBranch)
set_device_model!(template_ed, TapTransformer, StaticBranch)

set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_device_model!(template_ed, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_ed, PowerLoad, StaticPowerLoad)
set_device_model!(template_ed, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_ed, RenewableNonDispatch, FixedOutput)

set_device_model!(
        template_ed,
        DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless),  # HVDCTwoTerminalUnbounded # HVDCTwoTerminalLossless
    )


 problem_ed = DecisionModel(
    template_ed, 
    sys_ed; 
    optimizer = solver_ipopt,
    optimizer_solve_log_print = true, 
    name = "ED")

models = SimulationModels(;
        decision_models = [
            problem_uc,
            problem_ed,
        ],
    )

sequence = SimulationSequence(;
    models = models,
    feedforwards = Dict(
        "ED" => [
            SemiContinuousFeedforward(;
                component_type = ThermalStandard,
                source = OnVariable,
                affected_values = [ActivePowerVariable],
            ),
        ],
    ),
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(;
    name = "no_cache",
    steps = 2,
    models = models,
    sequence = sequence,
    simulation_folder = mktempdir(),
)

build_out = build!(sim)
@assert build_out == PSI.SimulationBuildStatus.BUILT

exports = Dict(
    "models" => [
        Dict(
            "name" => "UC",
            "store_all_variables" => true,
            "store_all_parameters" => true,
            "store_all_duals" => true,
            "store_all_aux_variables" => true,
        ),
        Dict(
            "name" => "ED",
            "store_all_variables" => true,
            "store_all_parameters" => true,
            "store_all_duals" => true,
            "store_all_aux_variables" => true,
        ),
    ],
    "path" => mktempdir(),
    "optimizer_stats" => true,
)
execute_out = execute!(sim; exports = exports, in_memory = true)
@assert execute_out == PSI.RunStatus.SUCCESSFULLY_FINALIZED

results = SimulationResults(sim);
uc_results = get_decision_problem_results(results, "UC")
ed_results = get_decision_problem_results(results, "ED")

# res_vars = read_variables(uc_results)
res_vars = read_variables(ed_results)

@show first(res_vars["FlowActivePowerFromToVariable__Line"])

@show first(res_vars["FlowActivePowerVariable__TwoTerminalGenericHVDCLine"])