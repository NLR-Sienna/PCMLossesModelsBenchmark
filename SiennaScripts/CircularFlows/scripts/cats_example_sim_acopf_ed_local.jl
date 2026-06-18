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
repo_path = joinpath(this_path, "..", "..", "..")
Pkg.activate(repo_path)
Pkg.instantiate()
# Uncomment this line to develop the local HSL_jll package (with MA57) for Ipopt's linear solver
#Pkg.develop(path = joinpath(repo_path, "lbt_HSL_jll.jl-2023.11.7", "HSL_jll.jl-2023.11.7"))
using HSL_jll

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
using Xpress
const PNM = PowerNetworkMatrices

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
xpress_milp = PSI.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.05)
ipopt_nlp  = JuMP.optimizer_with_attributes(() -> Ipopt.Optimizer(),
    "print_level" => 5,
    "hsllib" => HSL_jll.libhsl_path,
    "linear_solver" => "ma57"
)

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
    bounded = false,
    voltage_threshold = 275.0,
)
    ptdf = PTDF(sys)

    sim = build_uc_ed_simulation_with_acopf(
        sys, sys;
        uc_models    = build_cats_uc_models_hv(; voltage_threshold = voltage_threshold, bounded = bounded),
        ed_models    = build_cats_ed_models_acopf(; bounded = bounded),
        uc_optimizer = xpress_milp,
        ed_optimizer = ipopt_nlp,
        ptdf_uc      = ptdf,
        ptdf_ed      = ptdf,
        ignore_pf_ed = ignore_pf_ed,
        initial_time = DateTime("2019-01-01T12:00:00")
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
sys_a1 = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a1, 1.0)
scale_cats_loads!(sys_a1, 0.35)
# Running Unbounded with ignore_pf_ed=false leads to numerical issues in the AC power flow
C_a1, sim_a1 = detect_cats_circular_flows_sim_acopf(sys_a1; ignore_pf_ed = false, voltage_threshold = 275.0, bounded = true);
println("  Cycles found: $(length(C_a1))")
for (i, c) in enumerate(C_a1)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

sim_res_a = SimulationResults(sim_a1)
res_ed_a  = get_decision_problem_results(sim_res_a, "ED")
res_uc_a = get_decision_problem_results(sim_res_a, "UC")
uc_dual = read_realized_variable(res_uc_a, "CopperPlateBalanceConstraint__System")
stab_factors_a = read_realized_aux_variable(
    res_ed_a,
    "PowerFlowVoltageStabilityFactors__ACBus";
    table_format = TableFormat.WIDE,
)



flow_line_acopf = read_realized_variable(res_ed_a, "FlowActivePowerFromToVariable__Line")
flow_line_pf_aux = read_realized_aux_variable(res_ed_a, "PowerFlowBranchActivePowerFromTo__Line")
flow_hvdc_acopf = read_realized_variable(res_ed_a, "FlowReactivePowerFromToVariable__TwoTerminalGenericHVDCLine")
flow_hvdc_pf_aux = read_realized_aux_variable(res_ed_a, "PowerFlowBranchActivePowerFromTo__TwoTerminalGenericHVDCLine")


# ── Scenario B: no circular flow expected ──────

println()
println("=== Scenario B: no circular flow expected(ignore_pf_ed=false) ===")
sys_b = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_b, -1.0)
scale_cats_loads!(sys_b, 0.35)
C_b, sim_b = detect_cats_circular_flows_sim_acopf(sys_b; ignore_pf_ed = true, voltage_threshold = 275.0, bounded = true);
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end

sim_res_b = SimulationResults(sim_b)
res_ed_b  = get_decision_problem_results(sim_res_b, "ED")
res_uc_b = get_decision_problem_results(sim_res_b, "UC")
uc_dual_b = read_realized_variable(res_uc_b, "CopperPlateBalanceConstraint__System")
stab_factors_b = read_realized_aux_variable(
    res_ed_b,
    "PowerFlowVoltageStabilityFactors__ACBus";
    table_format = TableFormat.WIDE,
)

include(joinpath(this_path, "..", "print_utils.jl"))
print_stability_comparison(stab_factors_a, stab_factors_b; top_n = 20)

#=
┌───────────────────────┬───────────┬───────────────────────────┐
│ name                  │ available │ arc                       │
├───────────────────────┼───────────┼───────────────────────────┤
│ Newark_NRS_HVDC       │ true      │ Arc: bus-8335 -> bus-8814 │
│ Metcalf_SanJoseB_HVDC │ true      │ Arc: bus-1819 -> bus-1258 │
└───────────────────────┴───────────┴───────────────────────────┘
=#