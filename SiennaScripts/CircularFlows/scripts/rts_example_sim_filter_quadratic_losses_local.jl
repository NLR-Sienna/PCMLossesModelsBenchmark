# RTS cascaded UC (lossless PTDF, HiGHS) + ED (quadratic P=I²R losses, Ipopt)
# with low-voltage branch filtering to validate the filter approach.
#
# Lines and tap-transformers whose from-bus base voltage is at or below
# `voltage_threshold` kV are excluded from NetworkFlowConstraint.
# The FULL system PTDF is still used for the quadratic loss term.
#
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/rts_example_sim_filter_quadratic_losses_local.jl
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

function print_rts_branch_filter_stats(sys, voltage_threshold::Float64)
    n_lines    = length(collect(PSY.get_components(PSY.Line, sys)))
    n_hv_lines = count(
        x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold,
        PSY.get_components(PSY.Line, sys),
    )
    n_tap    = length(collect(PSY.get_components(PSY.TapTransformer, sys)))
    n_hv_tap = count(
        x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold,
        PSY.get_components(PSY.TapTransformer, sys),
    )
    println("  Branch filter (>$(voltage_threshold) kV from-bus): " *
            "Lines $(n_hv_lines)/$(n_lines), " *
            "TapTransformer $(n_hv_tap)/$(n_tap)")
end

"""
    detect_rts_circular_flows_sim_quad_filtered(
        sys;
        time_step = 1,
        voltage_threshold = 100.0
    ) -> Tuple{Vector{CircularFlow}, Simulation}

Build and execute a cascaded UC+ED simulation for the RTS system with filtered
branch flow constraints, then detect circular flows using PTDF-based graph
construction.

`ignore_pf_uc = true` and `ignore_pf_ed = true` are both hardcoded.  This is
required when branches are filtered: filtered branches have no AC PF aux variable
slots, so `build_graph_from_pf_aux_variables` would yield an incomplete graph.
The PTDF × bus-injection path is used instead.
"""
function detect_rts_circular_flows_sim_quad_filtered(
    sys;
    time_step::Int = 1,
    voltage_threshold::Float64 = 100.0,
)
    ptdf = PTDF(sys)

    print_rts_branch_filter_stats(sys, voltage_threshold)

    uc_models_filtered = build_rts_uc_models_hv(; voltage_threshold)
    ed_models_filtered = build_rts_ed_models_hv(; voltage_threshold)

    sim = build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(
        sys, sys;
        uc_models    = uc_models_filtered,
        ed_models    = ed_models_filtered,
        uc_optimizer = highs_milp,
        ed_optimizer = ipopt_nlp,
        ptdf_uc      = ptdf,
        ptdf_ed      = ptdf,
        ignore_pf_uc = true,
        ignore_pf_ed = true,
    )

    ed_container = sim.models.decision_models[2].internal.container

    execute!(sim; enable_progress_bar = false)

    optimize!(PSI.get_jump_model(ed_container))

    inj_expr = ed_container.expressions[
        PSY.InfrastructureSystems.Optimization.ExpressionKey{
            PSI.ActivePowerBalance, PSY.ACBus}("")
    ]
    inj_pu = JuMP.value.(inj_expr[:, 1]).data

    sim_res     = SimulationResults(sim)
    res_ed      = get_decision_problem_results(sim_res, "ED")
    res_vars_ed = read_hvdc_flow_variables(res_ed, read_realized_variable)

    G, bus_lookup = build_graph_from_ptdf(inj_pu, ptdf; base_power = PSY.get_base_power(sys))
    add_hvdc_edges!(G, sys, res_vars_ed, bus_lookup; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches), sim
end

println("=== RTS branch filter survey ===")
sys_survey = build_rts_system()
for thr in [0.0, 100.0, 138.0, 230.0, 345.0]
    n_lines = length(collect(PSY.get_components(PSY.Line, sys_survey)))
    n_hv    = count(x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > thr,
                    PSY.get_components(PSY.Line, sys_survey))
    n_tap   = length(collect(PSY.get_components(PSY.TapTransformer, sys_survey)))
    n_hv_t  = count(x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > thr,
                    PSY.get_components(PSY.TapTransformer, sys_survey))
    println("  threshold=$(thr) kV: Lines $(n_hv)/$(n_lines), TapTransformer $(n_hv_t)/$(n_tap)")
end
println()

println("=== Scenario A: with circular flow (positive RE cost_sign, filtered quadratic losses) ===")
sys_a = build_rts_system()
set_renewable_costs!(sys_a, 10.0)
C_a, sim_a = detect_rts_circular_flows_sim_quad_filtered(sys_a; voltage_threshold = 100.0)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

println("=== Scenario B: without circular flow (negative RE cost_sign, filtered quadratic losses) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0)
C_b, sim_b = detect_rts_circular_flows_sim_quad_filtered(sys_b; voltage_threshold = 100.0)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end
