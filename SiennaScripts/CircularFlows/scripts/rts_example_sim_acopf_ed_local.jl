# RTS cascaded UC (lossless PTDF + HV filter_function, HiGHS) +
# ED (ACPPowerModel all-lines, Ipopt).
# Two selectable modes via ignore_pf_ed:
#   true  (default) — circular flows from ACPPowerModel optimization variables directly
#   false           — post-solve AC power flow; circular flows from PF aux variables;
#                     voltage stability factors available from PSI results
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/rts_example_sim_acopf_ed_local.jl
using Pkg
this_path = @__DIR__
Pkg.activate(joinpath(this_path, "..", "..", ".."))
Pkg.instantiate()

using PowerSystems
using InfrastructureSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using PowerNetworkMatrices
const PNM = PowerNetworkMatrices
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
include(joinpath(this_path, "..", "..", "..", "Systems", "RTS", "build_rts.jl"))
include(joinpath(this_path, "..", "..", "build_models.jl"))
include(joinpath(this_path, "..", "..", "build_simulations.jl"))
include(joinpath(this_path, "..", "..", "utils.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

highs_milp = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
ipopt_nlp  = PSI.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)

"""
    detect_rts_circular_flows_sim_acopf(
        sys;
        time_step = 1,
        ignore_pf_ed = true,
    ) -> Vector{CircularFlow}

Build and execute a cascaded UC (lossless PTDF) + ED (ACPPowerModel) simulation
for RTS, then detect circular flows.

`ignore_pf_ed = true` (default): pure ACPPowerModel — no post-solve AC power flow.
  Circular flows are read from `FlowActivePowerFromToVariable__Line` / `__TapTransformer`
  optimization variables via `build_graph_from_acopf_variables`.

`ignore_pf_ed = false`: post-solve AC power flow runs after each ED solve.
  Circular flows are read from `PowerFlowBranchActivePowerFromTo__Line` PF aux variables
  via `build_graph_from_pf_aux_variables`. Voltage stability factors are available in
  `res_ed` (accessible after the function returns via `SimulationResults`).
"""
function detect_rts_circular_flows_sim_acopf(
    sys::PSY.System;
    time_step::Int = 1,
    ignore_pf_ed::Bool = true,
)
    ptdf = PTDF(sys)

    sim = build_uc_ed_simulation_with_acopf(
        sys, sys;
        uc_models    = build_rts_uc_models_hv(),
        ed_models    = build_rts_ed_models_acopf(),
        uc_optimizer = highs_milp,
        ed_optimizer = ipopt_nlp,
        ptdf_uc      = ptdf,
        ptdf_ed      = ptdf,
        ignore_pf_ed = ignore_pf_ed,
        initial_time = DateTime("2020-01-01T00:00:00")
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
    return find_circular_flows(G, bus_lookup, branches), sim
end

# ── Scenario A: circular flow expected (positive RE cost, HVDC enabled) ─────────

println("=== Scenario A / mode=variables (ignore_pf_ed=false): ACPPowerModel variables ===")
sys_a1 = build_rts_system()

# Exploring network reduction data
ptdf = PTDF(sys_a1)
PNM.populate_branch_maps_by_type!(ptdf.network_reduction_data)
ptdf.network_reduction_data.name_to_arc_map[Line]

set_renewable_costs!(sys_a1, 1.0)
C_a1, sim_a1 = detect_rts_circular_flows_sim_acopf(sys_a1; ignore_pf_ed = true);
println("  Cycles found: $(length(C_a1))")
for (i, c) in enumerate(C_a1)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

sim_a1_res = SimulationResults(sim_a1)
res_ed_a1  = get_decision_problem_results(sim_a1_res, "ED")
res_uc_a1 = get_decision_problem_results(sim_a1_res, "UC")
uc_dual = read_realized_variable(res_uc_a1, "CopperPlateBalanceConstraint__System")
line_flow_pf_aux = read_realized_aux_variable(res_ed_a1, "PowerFlowBranchActivePowerFromTo__Line")
line_flow_acopf = read_realized_variable(res_ed_a1, "FlowActivePowerFromToVariable__Line")

# ── Scenario B: no circular flow expected (HVDC disabled, negative RE cost) ──────

println("=== Scenario B: no circular flow (ignore_pf_ed=true) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -100.0)
C_b, sim_b = detect_rts_circular_flows_sim_acopf(sys_b; ignore_pf_ed = true);
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end


sim_b_res = SimulationResults(sim_b)
res_ed_b  = get_decision_problem_results(sim_b_res, "ED")
res_uc_b = get_decision_problem_results(sim_b_res, "UC")
uc_dual = read_realized_variable(res_uc_b, "CopperPlateBalanceConstraint__System")