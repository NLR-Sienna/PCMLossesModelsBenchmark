# Test: RTS filtered-line quadratic-losses simulation
#
# Validates:
#   1. build_rts_uc/ed_models_hv return DeviceModel objects with filter_function
#   2. filter_function passes all RTS lines at threshold=100.0 (all are ≥138 kV)
#   3. filter_function excludes 138 kV lines at threshold=230.0
#   4. Filtered simulation Scenario A (positive RE costs) finds ≥1 circular flow
#   5. Filtered simulation Scenario B (negative RE costs) finds 0 circular flows
#   6. Aggressive filter (threshold=230.0) still runs end-to-end without error
#
# CATS is NOT tested here — it takes too long to solve.
# threshold=100.0 passes all RTS branches (all ≥138 kV) → baseline for unfiltered behavior.
# threshold=230.0 filters the 138 kV lines → tests the filter path with actual exclusions.
#
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/test_rts_filter_quadratic_losses.jl
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
ipopt_nlp  = PSI.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

function run_test(name::String, cond::Bool)
    if cond
        println("  PASS: $name")
    else
        error("  FAIL: $name")
    end
end

# Internal helper — identical logic to detect_rts_circular_flows_sim_quad_filtered in Task 6
# but defined inline so this script is self-contained (including Task 6 would trigger its top-level runs).
function _run_rts_filtered_sim(sys; voltage_threshold::Float64 = 100.0, time_step::Int = 1)
    ptdf      = PTDF(sys)
    uc_models = build_rts_uc_models_hv(; voltage_threshold)
    ed_models = build_rts_ed_models_hv(; voltage_threshold)
    sim = build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(
        sys, sys;
        uc_models    = uc_models,
        ed_models    = ed_models,
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
    inj_pu      = JuMP.value.(inj_expr[:, 1]).data
    sim_res     = SimulationResults(sim)
    res_ed      = get_decision_problem_results(sim_res, "ED")
    res_vars_ed = read_hvdc_flow_variables(res_ed, read_realized_variable)
    G, bus_lookup = build_graph_from_ptdf(inj_pu, ptdf; base_power = PSY.get_base_power(sys))
    add_hvdc_edges!(G, sys, res_vars_ed, bus_lookup; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches)
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: factory function structure (no solver)
# ─────────────────────────────────────────────────────────────────────────────
println("=== Test 1: factory function structure ===")

uc_m = build_rts_uc_models_hv(; voltage_threshold = 100.0)
ed_m = build_rts_ed_models_hv(; voltage_threshold = 100.0)

run_test("UC models is a Dict",                 uc_m isa Dict)
run_test("UC Line → PSI.DeviceModel",           uc_m[Line] isa PSI.DeviceModel)
run_test("UC TapTransformer → PSI.DeviceModel", uc_m[TapTransformer] isa PSI.DeviceModel)
run_test("UC Line has filter_function",         haskey(uc_m[Line].attributes, "filter_function"))
run_test("UC TapTransformer has filter_function", haskey(uc_m[TapTransformer].attributes, "filter_function"))
run_test("ED models is a Dict",                 ed_m isa Dict)
run_test("ED Line → PSI.DeviceModel",           ed_m[Line] isa PSI.DeviceModel)
run_test("ED TapTransformer → PSI.DeviceModel", ed_m[TapTransformer] isa PSI.DeviceModel)
println()

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: filter_function logic on real RTS branches (no solver)
# ─────────────────────────────────────────────────────────────────────────────
println("=== Test 2: filter_function behavior on real RTS branches ===")

sys_filter_test = build_rts_system()

filter_fn_100 = uc_m[Line].attributes["filter_function"]
all_pass_100 = all(filter_fn_100(br) for br in PSY.get_components(PSY.Line, sys_filter_test))
run_test("All RTS lines pass threshold=100.0 (all are ≥138 kV)", all_pass_100)

uc_m_230    = build_rts_uc_models_hv(; voltage_threshold = 230.0)
filter_fn_230 = uc_m_230[Line].attributes["filter_function"]
some_fail_230 = any(!filter_fn_230(br) for br in PSY.get_components(PSY.Line, sys_filter_test))
run_test("Some RTS lines fail threshold=230.0 (138 kV lines excluded)", some_fail_230)
println()

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Scenario A — positive RE costs, threshold=100.0 → ≥1 cycle expected
# ─────────────────────────────────────────────────────────────────────────────
println("=== Test 3: Scenario A — circular flow expected (threshold=100.0) ===")

sys_a = build_rts_system()
set_renewable_costs!(sys_a, 10.0)
C_a = _run_rts_filtered_sim(sys_a; voltage_threshold = 100.0)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("    Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end
run_test("Scenario A finds ≥1 circular flow", length(C_a) >= 1)
println()

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: Scenario B — negative RE costs, threshold=100.0 → 0 cycles expected
# ─────────────────────────────────────────────────────────────────────────────
println("=== Test 4: Scenario B — no circular flow expected (threshold=100.0) ===")

sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0)
C_b = _run_rts_filtered_sim(sys_b; voltage_threshold = 100.0)
println("  Cycles found: $(length(C_b))")
run_test("Scenario B finds 0 circular flows", length(C_b) == 0)
println()

# ─────────────────────────────────────────────────────────────────────────────
# Test 5: Aggressive filter (threshold=230.0) — simulation completes without error
# ─────────────────────────────────────────────────────────────────────────────
println("=== Test 5: Scenario A with aggressive filter (threshold=230.0) ===")

sys_c = build_rts_system()
set_renewable_costs!(sys_c, 10.0)
C_c = _run_rts_filtered_sim(sys_c; voltage_threshold = 230.0)
println("  Cycles found: $(length(C_c))")
run_test("Aggressive filter (230.0 kV) completes without error", true)
println()

println("=== All tests passed! ===")
