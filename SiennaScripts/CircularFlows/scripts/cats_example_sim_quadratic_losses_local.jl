# CATS cascaded UC (lossless PTDF, HiGHS) + ED (quadratic P=I²R losses, Ipopt)
# using PSI Simulation with SemiContinuousFeedforward — no manual binary fixing.
# Binaries are propagated from UC → ED via feedforward; ED is a pure NLP.
# No Gurobi required.
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/cats_example_sim_quadratic_losses_local.jl
using Pkg
this_path = @__DIR__
Pkg.activate(joinpath(this_path, "..", "..", ".."))
Pkg.instantiate()

using PowerSystems
using InfrastructureSystems
using PowerSimulations
using HydroPowerSimulations
using PowerFlows
using PowerNetworkMatrices
using Graphs
using SimpleWeightedGraphs
using JuMP
using HiGHS
using Ipopt
using Dates
using DataFrames
using Logging

include(joinpath(this_path, "..", "mapped_indices.jl"))
include(joinpath(this_path, "..", "circular_flows.jl"))
include(joinpath(this_path, "..", "..", "..", "SiennaScripts", "add_hvdc.jl"))
include(joinpath(this_path, "..", "..", "..", "Systems", "CATS", "build_cats.jl"))
include(joinpath(this_path, "..", "..", "build_models.jl"))
include(joinpath(this_path, "..", "..", "build_simulations.jl"))
include(joinpath(this_path, "..", "..", "utils.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

cats_json = joinpath(this_path, "..", "..", "..", "Systems", "CATS", "CATS_saved_reduced_sys.json")

highs_milp = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
ipopt_nlp  = PSI.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)

"""
    detect_cats_circular_flows_sim_quad(
        sys;
        time_step::Int = 1
    ) -> Vector{CircularFlow}

Build and execute a cascaded UC+ED simulation for the CATS system using PTDF-based
graph construction, then detect circular flows in the ED solution.

`ignore_pf_uc = true` and `ignore_pf_ed = true` are both hardcoded: neither stage
runs the post-solve AC power flow, so graph construction always uses PTDF × bus-injection
estimates rather than actual AC branch flows.

# Arguments
- `sys`: `PSY.System` to optimise (modified in-place by HVDC additions upstream).
- `time_step`: which time-step of the simulation to analyse (default: `1`).

# Returns
- `Vector{CircularFlow}` — detected closed loops and their branch flows.
"""
function detect_cats_circular_flows_sim_quad(sys; time_step::Int = 1)
    ptdf = PTDF(sys)

    # filt_func(x) = get_base_voltage(get_from(get_arc(x))) > 325.0
    # Filter functions are passed as attributes to set_device_model
    # example: set_device_model!(template_uc, DeviceModel(Line, StaticBranch; attributes = Dict("filter_function" => x -> get_base_voltage(get_from(get_arc(x))) > 100)))
    # Only include high-voltage branches in the graph, so we don't find trivial cycles with low-voltage distribution lines.


    # Build UC (lossless, HiGHS) + ED (quadratic losses, Ipopt) simulation.
    # SemiContinuousFeedforward propagates UC on/off to ED — no manual binary fixing.
    # ignore_pf_uc and ignore_pf_ed are both true: neither stage runs post-solve AC power
    # flow, so flows are inferred via PTDF × injections rather than AC branch flows.
    sim = build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(
        sys, sys;
        uc_models    = CATS_UC_MODELS,
        ed_models    = CATS_ED_MODELS,
        uc_optimizer = highs_milp,
        ed_optimizer = ipopt_nlp,
        ptdf_uc      = ptdf,
        ptdf_ed      = ptdf,
        ignore_pf_uc = true,
        ignore_pf_ed = true,
    )

    ed_container = sim.models.decision_models[2].internal.container

    execute!(sim; enable_progress_bar = false)

    # PSI clears the JuMP solution after execute!; re-optimize to restore it.
    optimize!(PSI.get_jump_model(ed_container))

    # Get bus injections from the ED container (per-unit, used for PTDF flows).
    inj_expr = ed_container.expressions[
        PSY.InfrastructureSystems.Optimization.ExpressionKey{
            PSI.ActivePowerBalance, PSY.ACBus}("")
    ]
    inj_pu = JuMP.value.(inj_expr[:, 1]).data

    # Get HVDC flows from simulation results for graph construction.
    # Uses read_realized_variable (not read_variable) for SimulationProblemResults.
    sim_res = SimulationResults(sim)
    res_ed  = get_decision_problem_results(sim_res, "ED")
    res_vars_ed = read_hvdc_flow_variables(res_ed, read_realized_variable)

    G, bus_lookup = build_graph_from_ptdf(inj_pu, ptdf; base_power = PSY.get_base_power(sys))
    add_hvdc_edges!(G, sys, res_vars_ed, bus_lookup; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches)
end

println("=== Scenario A: with circular flow (positive RE cost_sign, quadratic losses in ED) ===")
sys_a = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a, 1.0)
C_a = detect_cats_circular_flows_sim_quad(sys_a)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

println("=== Scenario B: without circular flow (HVDC disabled) ===")
sys_b = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_b, -1.0)
for hvdc in PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys_b)
    PSY.set_available!(hvdc, false)
end
C_b = detect_cats_circular_flows_sim_quad(sys_b)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

### Double check if there is congestion when having circular flows in positive renewable cost scenario
### Avoid re-running optimize to get the PTDF injections.
### Filter lines with losses so don't include every line in the model ###
### Try to run CATS with full OPF in ED ###