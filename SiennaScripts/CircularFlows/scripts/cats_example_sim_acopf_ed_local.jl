# CATS cascaded UC (lossless PTDF + HV filter_function, HiGHS) +
# ED (ACPPowerModel all-lines, Ipopt).
# Two selectable modes via ignore_pf_ed:
#   true  (default) — circular flows from ACPPowerModel optimization variables directly
#   false           — post-solve AC power flow; circular flows from PF aux variables;
#                     voltage stability factors available from PSI results
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/cats_example_sim_acopf_ed_local.jl
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
    detect_cats_circular_flows_sim_acopf(
        sys;
        time_step = 1,
        ignore_pf_ed = true,
    ) -> Vector{CircularFlow}

Build and execute a cascaded UC (lossless PTDF) + ED (ACPPowerModel) simulation
for CATS, then detect circular flows.

`ignore_pf_ed = true` (default): pure ACPPowerModel — no post-solve AC power flow.
  Circular flows are read from `FlowActivePowerFromToVariable__Line` / `__Transformer2W`
  optimization variables via `build_graph_from_acopf_variables`.

`ignore_pf_ed = false`: post-solve AC power flow runs after each ED solve.
  Circular flows are read from `PowerFlowBranchActivePowerFromTo__Line` PF aux variables
  via `build_graph_from_pf_aux_variables`. Voltage stability factors are available in
  `res_ed` (accessible after the function returns via `SimulationResults`).
"""
function detect_cats_circular_flows_sim_acopf(
    sys::PSY.System;
    time_step::Int = 1,
    ignore_pf_ed::Bool = true,
)
    ptdf = PTDF(sys)

    sim = build_uc_ed_simulation_with_acopf(
        sys, sys;
        uc_models    = build_cats_uc_models_hv(),
        ed_models    = build_cats_ed_models_acopf(),
        uc_optimizer = highs_milp,
        ed_optimizer = ipopt_nlp,
        ptdf_uc      = ptdf,
        ptdf_ed      = ptdf,
        ignore_pf_ed = ignore_pf_ed,
    )

    execute!(sim; enable_progress_bar = false)

    sim_res = SimulationResults(sim)
    res_ed  = get_decision_problem_results(sim_res, "ED")

    if ignore_pf_ed
        G, bus_lookup = build_graph_from_acopf_variables(res_ed, sys; time_step = time_step)
    else
        G, bus_lookup = build_graph_from_pf_aux_variables(res_ed, sys; time_step = time_step)
    end

    res_vars_ed = read_hvdc_flow_variables(res_ed, read_realized_variable)
    add_hvdc_edges!(G, sys, res_vars_ed, bus_lookup; time_step = time_step)

    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches)
end

# ── Scenario A: circular flow expected (positive RE cost, HVDC enabled) ─────────

println("=== Scenario A / mode=variables (ignore_pf_ed=true): ACPPowerModel variables ===")
sys_a1 = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a1, 1.0)
C_a1 = detect_cats_circular_flows_sim_acopf(sys_a1; ignore_pf_ed = true)
println("  Cycles found: $(length(C_a1))")
for (i, c) in enumerate(C_a1)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

println()
println("=== Scenario A / mode=pf (ignore_pf_ed=false): post-solve AC PF ===")
println("    (voltage stability factors available in res_ed after this call)")
sys_a2 = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a2, 1.0)
C_a2 = detect_cats_circular_flows_sim_acopf(sys_a2; ignore_pf_ed = false)
println("  Cycles found: $(length(C_a2))")
for (i, c) in enumerate(C_a2)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

# ── Scenario B: no circular flow expected (HVDC disabled, negative RE cost) ──────

println()
println("=== Scenario B: no circular flow (HVDC disabled, ignore_pf_ed=true) ===")
sys_b = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_b, -1.0)
for hvdc in PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys_b)
    PSY.set_available!(hvdc, false)
end
C_b = detect_cats_circular_flows_sim_acopf(sys_b; ignore_pf_ed = true)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end
