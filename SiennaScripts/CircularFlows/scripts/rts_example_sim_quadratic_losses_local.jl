# RTS cascaded UC (lossless PTDF, HiGHS) + ED (quadratic P=I²R losses, Ipopt)
# using PSI Simulation with SemiContinuousFeedforward — no manual binary fixing.
# Binaries are propagated from UC → ED via feedforward; ED is a pure NLP.
# No Gurobi required.
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/rts_example_sim_quadratic_losses_local.jl
using Pkg
this_path = @__DIR__
Pkg.activate(joinpath(this_path, "..", "..", ".."))
Pkg.instantiate()

using PowerSystemCaseBuilder
using PowerSystems
using InfrastructureSystems
using PowerSimulations
using HydroPowerSimulations
using PowerFlows
using PowerNetworkMatrices
using Graphs
using SimpleWeightedGraphs
using HiGHS
using Ipopt
using Dates
using JuMP
using DataFrames
using Logging

include(joinpath(this_path, "..", "mapped_indices.jl"))
include(joinpath(this_path, "..", "circular_flows.jl"))
include(joinpath(this_path, "..", "print_utils.jl"))
include(joinpath(this_path, "..", "..", "..", "Systems", "RTS", "build_rts.jl"))
include(joinpath(this_path, "..", "..", "build_models.jl"))
include(joinpath(this_path, "..", "..", "build_simulations.jl"))
include(joinpath(this_path, "..", "..", "utils.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

highs_milp = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
ipopt_nlp  = PSI.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)

# RTS-specific UC/ED device models.
# RTS uses TapTransformer (not Transformer2W) and no SynchronousCondenser.
const RTS_UC_MODELS = Dict(
    Line                        => StaticBranchBounds,
    TapTransformer              => StaticBranchBounds,
    ThermalStandard             => ThermalBasicUnitCommitment,
    PowerLoad                   => StaticPowerLoad,
    RenewableDispatch           => RenewableFullDispatch,
    HydroDispatch               => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine  => HVDCTwoTerminalLossless,
)

const RTS_ED_MODELS = Dict(
    Line                        => StaticBranchBounds,
    TapTransformer              => StaticBranchBounds,
    ThermalStandard             => ThermalBasicDispatch,
    PowerLoad                   => StaticPowerLoad,
    RenewableDispatch           => RenewableFullDispatch,
    HydroDispatch               => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine  => HVDCTwoTerminalLossless,
)

"""
    detect_rts_circular_flows_sim_quad(
        sys;
        time_step::Int = 1,
        use_pf::Bool = true
    ) -> Tuple{Vector{CircularFlow}, Simulation}

Build and execute a cascaded UC+ED simulation for the RTS system, then detect
circular flows in the ED solution.

`ignore_pf_uc = true` is hardcoded: UC never runs the post-solve AC power flow.
Only the ED stage runs AC power flow post-solve, so voltage stability factors are
always available regardless of the `use_pf` setting.

# Arguments
- `sys`: `PSY.System` to optimise.
- `time_step`: which time-step of the simulation to analyse (default: `1`).
- `use_pf`: when `true` (default), builds the flow graph from actual AC power
  flow branch flows stored as PSI auxiliary variables; when `false`, falls back
  to PTDF × bus-injection estimates. Both paths still produce voltage stability
  factors because the ED always runs AC power flow post-solve.

# Returns
- `(circular_flows::Vector{CircularFlow}, sim::Simulation)`
"""
function detect_rts_circular_flows_sim_quad(sys; time_step::Int = 1, use_pf::Bool = true)
    ptdf = PTDF(sys)

    # ED always runs AC power flow post-solve so that voltage stability factors
    # are available regardless of which graph-building method is chosen.
    sim = build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(
        sys, sys;
        uc_models    = RTS_UC_MODELS,
        ed_models    = RTS_ED_MODELS,
        uc_optimizer = highs_milp,
        ed_optimizer = ipopt_nlp,
        ptdf_uc      = ptdf,
        ptdf_ed      = ptdf,
        ignore_pf_uc = true,
        ignore_pf_ed = false,
    )

    execute!(sim; enable_progress_bar = false)

    # Retrieve ED results for graph construction and stability factor extraction.
    sim_res = SimulationResults(sim)
    res_ed  = get_decision_problem_results(sim_res, "ED")

    # Get HVDC flows from simulation results for graph construction.
    res_vars_ed = read_hvdc_flow_variables(res_ed, read_realized_variable)

    if use_pf
        # Build graph from actual AC power flow branch flows (default).
        G, bus_lookup = build_graph_from_pf_aux_variables(res_ed, sys; time_step = time_step)
    else
        # Fall back to PTDF × injection flows. PSI clears the JuMP solution after
        # execute!, so re-optimize to restore variable values before reading them.
        ed_container = sim.models.decision_models[2].internal.container
        optimize!(PSI.get_jump_model(ed_container))
        inj_expr = ed_container.expressions[
            PSY.InfrastructureSystems.Optimization.ExpressionKey{
                PSI.ActivePowerBalance, PSY.ACBus}("")
        ]
        inj_pu = JuMP.value.(inj_expr[:, 1]).data
        G, bus_lookup = build_graph_from_ptdf(inj_pu, ptdf; base_power = PSY.get_base_power(sys))
    end

    add_hvdc_edges!(G, sys, res_vars_ed, bus_lookup; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches), sim
end

println("=== Scenario A: with circular flow (positive RE cost_sign, quadratic losses in ED) ===")
sys_a = build_rts_system()
set_renewable_costs!(sys_a, 10.0)
C_a, sim_a = detect_rts_circular_flows_sim_quad(sys_a; use_pf = false)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

sim_res_a = SimulationResults(sim_a)
res_ed_a  = get_decision_problem_results(sim_res_a, "ED")
stab_factors_a = read_realized_aux_variable(
    res_ed_a,
    "PowerFlowVoltageStabilityFactors__ACBus";
    table_format = TableFormat.WIDE,
)

println("=== Scenario B: without circular flow (negative RE cost_sign, quadratic losses in ED) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0)
C_b, sim_b = detect_rts_circular_flows_sim_quad(sys_b; use_pf = false)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

sim_res_b = SimulationResults(sim_b)
res_ed_b  = get_decision_problem_results(sim_res_b, "ED")
stab_factors_b = read_realized_aux_variable(
    res_ed_b,
    "PowerFlowVoltageStabilityFactors__ACBus";
    table_format = TableFormat.WIDE,
)

print_stability_comparison(stab_factors_a, stab_factors_b; top_n = 42)