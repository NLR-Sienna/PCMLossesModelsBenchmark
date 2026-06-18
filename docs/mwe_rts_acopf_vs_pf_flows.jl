# =============================================================================
# Minimum Working Example: ACOPF branch flows vs post-solve AC power-flow flows
# =============================================================================
#
# Purpose
# -------
# Build a small two-stage UC (PTDF) + ED (full AC OPF) simulation for the
# modified RTS-GMLC system (which includes a TwoTerminalGenericHVDCLine), then
# compare, branch by branch, two quantities that should be (nearly) identical:
#
#   1. FlowActivePowerFromToVariable__Line
#        -> the from->to active power flow taken DIRECTLY from the ACPPowerModel
#           (full AC OPF) optimization variables.
#
#   2. PowerFlowBranchActivePowerFromTo__Line
#        -> the from->to active power flow reported by the PowerFlows.jl
#           post-solve AC power flow (the auxiliary variable produced when the
#           ED NetworkModel is given `power_flow_evaluation = ACPowerFlow(...)`).
#
# Both are evaluated on the SAME dispatch solution (same bus injections coming
# out of the ED solve), so they should agree to within solver tolerance. They
# do not — the differences are large — which points at a bug in PowerFlows.jl
# (or in how PowerSimulations hands the ED solution to it). This script
# isolates and prints that discrepancy.
#
# This file is deliberately self-contained: it depends ONLY on PowerSystems,
# PowerSimulations, PowerSystemCaseBuilder, PowerFlows and the public solver /
# data packages. It does NOT use any helper from the repository it ships with,
# so the PowerFlows developer can run it directly.
#
# Circular-flow detection is intentionally OUT OF SCOPE here. The single thing
# this MWE highlights is the per-branch difference between the two flow values.
#
# How to run
# ----------
#   julia docs/mwe_rts_acopf_vs_pf_flows.jl
#
# Package versions this was reproduced with (Sienna registry):
#   PowerSystems          5.10.0
#   PowerSimulations      0.35.0
#   PowerFlows            0.18.0
#   PowerNetworkMatrices  0.21.1
#   PowerSystemCaseBuilder 2.2.1
#   HydroPowerSimulations 0.16.0
#   InfrastructureSystems 3.6.0
# =============================================================================

using Pkg

# Spin up a throwaway environment with just the packages this MWE needs, so it
# can be run anywhere without dragging along the host repo's Project.toml.
# Comment this block out and `using` your own environment if you prefer to pin
# exact versions via a Manifest you already have.
Pkg.activate(mktempdir())
Pkg.add([
    "PowerSystems",
    "PowerSimulations",
    "PowerSystemCaseBuilder",
    "PowerFlows",
    "PowerNetworkMatrices",
    "HydroPowerSimulations",
    "HiGHS",
    "Ipopt",
    "DataFrames",
])

using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using HiGHS
using Ipopt
using DataFrames
using Dates
using Logging

const PSY = PowerSystems
const PSI = PowerSimulations

highs_milp = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
ipopt_nlp  = PSI.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# -----------------------------------------------------------------------------
# 1. Build the modified RTS-GMLC system (with HVDC).
#
# This mirrors a typical RTS setup: hourly single-time-series, one HVDC line
# made lossless with symmetric ±100 MW P/Q limits, and one hydro unit disabled.
# -----------------------------------------------------------------------------
function build_rts_system_mwe()
    sys = build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization = true)
    transform_single_time_series!(sys, Hour(1), Hour(1))

    hvdc = only(get_components(TwoTerminalGenericHVDCLine, sys))
    set_loss!(hvdc, LinearCurve(0.0))
    set_reactive_power_limits_from!(hvdc, (min = -100, max = 100))
    set_reactive_power_limits_to!(hvdc, (min = -100, max = 100))
    set_active_power_limits_from!(hvdc, (min = -100, max = 100))
    set_active_power_limits_to!(hvdc, (min = -100, max = 100))

    set_available!(get_component(HydroDispatch, sys, "201_HYDRO_4"), false)
    return sys
end

# -----------------------------------------------------------------------------
# 2. UC template: lossless PTDF copper-plate. Plain StaticBranchBounds on the
#    branches (no LV filtering — that only matters for circular-flow work, which
#    is out of scope here). HVDC modeled lossless.
# -----------------------------------------------------------------------------
function make_uc_template(ptdf)
    template = ProblemTemplate(
        NetworkModel(
            PTDFPowerModel;
            PTDF_matrix = ptdf,
            use_slacks = true,
            duals = [CopperPlateBalanceConstraint],
        ),
    )
    set_device_model!(template, Line, StaticBranchBounds)
    set_device_model!(template, TapTransformer, StaticBranchBounds)
    set_device_model!(template, ThermalStandard, ThermalBasicUnitCommitment)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless)
    return template
end

# -----------------------------------------------------------------------------
# 3. ED template: full AC OPF (ACPPowerModel). All branches are
#    StaticBranchUnbounded so the AC OPF sees the complete network.
#
#    THE KEY LINE for this MWE:
#       power_flow_evaluation = PowerFlows.ACPowerFlow(...)
#    This asks PowerSimulations to run a PowerFlows.jl AC power flow on the ED
#    solution AFTER the solve, and store the resulting branch flows as the
#    auxiliary variable PowerFlowBranchActivePowerFromTo__Line. That aux
#    variable is what we compare against the ACOPF optimization variable
#    FlowActivePowerFromToVariable__Line.
# -----------------------------------------------------------------------------
function make_ed_template()
    template = ProblemTemplate(
        NetworkModel(
            ACPPowerModel;
            use_slacks = true,
            power_flow_evaluation = PowerFlows.ACPowerFlow(;
                calculate_loss_factors = true,
                calculate_voltage_stability_factors = true,
            ),
        ),
    )
    set_device_model!(template, Line, StaticBranchUnbounded)
    set_device_model!(template, TapTransformer, StaticBranchUnbounded)
    set_device_model!(template, ThermalStandard, ThermalBasicDispatch)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, TwoTerminalGenericHVDCLine, HVDCTwoTerminalLossless)
    return template
end

# -----------------------------------------------------------------------------
# 4. Assemble and run the UC -> ED simulation (1 step).
# -----------------------------------------------------------------------------
function build_and_run_sim(sys)
    ptdf = PTDF(sys)

    uc_model = DecisionModel(
        make_uc_template(ptdf), sys;
        optimizer = highs_milp, name = "UC", store_variable_names = true,
    )
    ed_model = DecisionModel(
        make_ed_template(), sys;
        optimizer = ipopt_nlp, name = "ED", store_variable_names = true,
    )

    models = SimulationModels(; decision_models = [uc_model, ed_model])

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
        name = "rts_acopf_vs_pf_mwe",
        steps = 1,
        models = models,
        sequence = sequence,
        simulation_folder = mktempdir(),
        initial_time = DateTime("2020-01-01T00:00:00"),
    )

    build!(sim; console_level = Logging.Error)
    execute!(sim; enable_progress_bar = false)
    return sim
end

# -----------------------------------------------------------------------------
# 5. Compare the two flow quantities, branch by branch, at one time step.
# -----------------------------------------------------------------------------
function compare_flows(sim; time_step::Int = 1)
    res = get_decision_problem_results(SimulationResults(sim), "ED")

    # ACOPF optimization variable (from->to active power) for Lines.
    acopf = read_realized_variable(res, "FlowActivePowerFromToVariable__Line")

    # PowerFlows.jl post-solve AC power-flow result (from->to active power) for Lines.
    pf_aux = read_realized_aux_variable(res, "PowerFlowBranchActivePowerFromTo__Line")

    # Both come back in long format with columns (DateTime, name, value).
    # Note the two tables can even have a DIFFERENT number of branch rows — the
    # PF aux table typically lists more branches than the ACOPF flow variable —
    # which is itself worth flagging. Pick one time step and inner-join on name.
    timestamps = sort(unique(acopf.DateTime))
    ts = timestamps[time_step]
    println("  comparing at timestamp ", ts,
            "  (ACOPF branches=", count(==(ts), acopf.DateTime),
            ", PF-aux branches=", count(==(ts), pf_aux.DateTime), ")")

    a = acopf[acopf.DateTime .== ts, [:name, :value]]
    p = pf_aux[pf_aux.DateTime .== ts, [:name, :value]]

    df = innerjoin(a, p, on = :name, makeunique = true)
    rename!(df, :name => :line, :value => :acopf_var_MW, :value_1 => :pf_aux_MW)
    df.abs_diff_MW = abs.(df.acopf_var_MW .- df.pf_aux_MW)

    sort!(df, :abs_diff_MW; rev = true)
    return df
end

# -----------------------------------------------------------------------------
# Run it.
# -----------------------------------------------------------------------------
println("Building RTS system (with HVDC)...")
sys = build_rts_system_mwe()

println("Building and running UC (PTDF) -> ED (AC OPF) simulation...")
sim = build_and_run_sim(sys)

println("\nComparing ACOPF flow variable vs PowerFlows.jl post-solve AC PF (Lines, t=1):")
df = compare_flows(sim; time_step = 1)

show(df, allrows = true, allcols = true)
println()
println("\nSummary:")
println("  branches compared : ", nrow(df))
println("  max |diff| (MW)   : ", round(maximum(df.abs_diff_MW); digits = 4))
println("  mean |diff| (MW)  : ", round(sum(df.abs_diff_MW) / nrow(df); digits = 4))
println("\nThe two columns above are computed on the SAME ED dispatch solution and")
println("should match to solver tolerance. The large differences are the issue to")
println("investigate in PowerFlows.jl.")
