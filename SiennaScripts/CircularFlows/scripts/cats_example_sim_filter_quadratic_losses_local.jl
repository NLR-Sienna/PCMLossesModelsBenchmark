# CATS cascaded UC (lossless PTDF, HiGHS) + ED (quadratic P=I²R losses, Ipopt)
# with low-voltage branch filtering to reduce model size.
#
# Lines and transformers whose from-bus base voltage is at or below
# `voltage_threshold` kV are excluded from NetworkFlowConstraint (thermal-limit
# rows).  The FULL system PTDF is still used for the quadratic loss term so that
# losses on ALL branches — including filtered LV lines — are correctly captured.
#
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/cats_example_sim_filter_quadratic_losses_local.jl
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

function print_branch_filter_stats(sys, voltage_threshold::Float64)
    n_lines    = length(collect(PSY.get_components(PSY.Line, sys)))
    n_hv_lines = count(
        x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold,
        PSY.get_components(PSY.Line, sys),
    )
    n_t2w    = length(collect(PSY.get_components(PSY.Transformer2W, sys)))
    n_hv_t2w = count(
        x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold,
        PSY.get_components(PSY.Transformer2W, sys),
    )
    println("  Branch filter (>$(voltage_threshold) kV from-bus): " *
            "Lines $(n_hv_lines)/$(n_lines), " *
            "Transformer2W $(n_hv_t2w)/$(n_t2w)")
end

"""
    detect_cats_circular_flows_sim_quad_filtered(
        sys;
        time_step = 1,
        voltage_threshold = 100.0
    ) -> Vector{CircularFlow}

Build and execute a cascaded UC+ED simulation for the CATS system with filtered
branch flow constraints, then detect circular flows in the ED solution.

Lines and transformers whose from-bus base voltage is ≤ `voltage_threshold` kV
are excluded from `NetworkFlowConstraint` to reduce MILP/NLP size.  The full
system PTDF is used for both the network model and the quadratic loss term so
that losses on all branches are captured correctly regardless of the filter.

`ignore_pf_uc = true` and `ignore_pf_ed = true` are both hardcoded: neither
stage runs a post-solve AC power flow.  This is required when branches are
filtered because filtered branches have no AC PF aux variable slots, making
`build_graph_from_pf_aux_variables` incomplete.  Graph construction uses
PTDF × bus-injection estimates instead.
"""
function detect_cats_circular_flows_sim_quad_filtered(
    sys;
    time_step::Int = 1,
    voltage_threshold::Float64 = 100.0,
)
    ptdf = PTDF(sys)

    print_branch_filter_stats(sys, voltage_threshold)

    uc_models_filtered = build_cats_uc_models_hv(; voltage_threshold)
    ed_models_filtered = build_cats_ed_models_hv(; voltage_threshold)

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
    return find_circular_flows(G, bus_lookup, branches)
end

println("=== CATS branch filter survey ===")
sys_survey = build_cats_system(cats_json)
for thr in [0.0, 69.0, 100.0, 230.0, 345.0]
    n_lines = length(collect(PSY.get_components(PSY.Line, sys_survey)))
    n_hv    = count(x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > thr,
                    PSY.get_components(PSY.Line, sys_survey))
    n_t2w   = length(collect(PSY.get_components(PSY.Transformer2W, sys_survey)))
    n_hv_t2 = count(x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > thr,
                    PSY.get_components(PSY.Transformer2W, sys_survey))
    println("  threshold=$(thr) kV: Lines $(n_hv)/$(n_lines), Transformer2W $(n_hv_t2)/$(n_t2w)")
end
println()

println("=== Scenario A: with circular flow (positive RE cost_sign, filtered quadratic losses) ===")
sys_a = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a, 1.0)
C_a = detect_cats_circular_flows_sim_quad_filtered(sys_a; voltage_threshold = 345.0)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

println("=== Scenario B: without circular flow (HVDC disabled, filtered quadratic losses) ===")
sys_b = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_b, -1.0)
for hvdc in PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys_b)
    PSY.set_available!(hvdc, false)
end
C_b = detect_cats_circular_flows_sim_quad_filtered(sys_b; voltage_threshold = 100.0)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end
