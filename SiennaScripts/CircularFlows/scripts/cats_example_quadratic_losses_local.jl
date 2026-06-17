# CATS cascaded UC (lossless PTDF, HiGHS) + ED (quadratic P=I²R losses, Ipopt).
# Binaries fixed by feedforward → ED is a pure NLP, no Gurobi required.
# Run from the repo root:
#   julia SiennaScripts/CircularFlows/scripts/cats_example_quadratic_losses_local.jl
using Pkg
this_path = @__DIR__
Pkg.activate(joinpath(this_path, "..", "..", ".."))
Pkg.instantiate()

using PowerSystems
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
include(joinpath(this_path, "..", "..", "..", "Systems", "CATS", "build_cats.jl"))
include(joinpath(this_path, "..", "..", "build_models.jl"))
include(joinpath(this_path, "..", "..", "utils.jl"))
include(joinpath(this_path, "..", "..", "..", "SiennaScripts", "add_hvdc.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

cats_json = joinpath(this_path, "..", "..", "..", "Systems", "CATS", "CATS_saved_reduced_sys.json")

highs_milp = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
ipopt_nlp  = PSI.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)

function get_ptdf_bus_injections(model)
    container = PSI.get_optimization_container(model)
    optimize!(PSI.get_jump_model(container))
    inj_expr = container.expressions[
        PSY.InfrastructureSystems.Optimization.ExpressionKey{
            PSI.ActivePowerBalance, PSY.ACBus}("")
    ]
    return JuMP.value.(inj_expr[:, 1]).data
end

function detect_cats_circular_flows_quad(sys; time_step::Int = 1)
    ptdf = PTDF(sys)

    # Step 1: lossless UC to get commitment decisions
    uc_model = make_ptdf_model_without_losses(
        sys;
        device_models = CATS_UC_MODELS,
        optimizer     = highs_milp,
        ptdf          = ptdf,
        name          = "UC",
        ignore_pf     = true,
    )
    PSI.build!(uc_model, output_dir = mktempdir())
    PSI.solve!(uc_model)

    # Step 2: ED with quadratic losses; fix UC binary decisions
    ed_model = make_ptdf_model_without_losses(
        sys;
        device_models = CATS_ED_MODELS,
        optimizer     = ipopt_nlp,
        ptdf          = ptdf,
        name          = "ED",
    )
    PSI.build!(ed_model, output_dir = mktempdir())

    uc_sol = Dict(
        name(v) => value(v)
        for v in all_variables(PSI.get_jump_model(PSI.get_optimization_container(uc_model)))
    )
    for v in all_variables(PSI.get_jump_model(PSI.get_optimization_container(ed_model)))
        if is_binary(v)
            val = get(uc_sol, name(v), 0.0)
            unset_binary(v)
            set_lower_bound(v, val)
            set_upper_bound(v, val)
        end
    end

    # Step 3: add quadratic losses to ED (flat voltage, no prior AC solve needed)
    update_copperplate_quadratic_loss_approximation_no_voltage!(ed_model, sys, ptdf)
    PSI.solve!(ed_model)

    results_ed = PSI.OptimizationProblemResults(ed_model)
    res_vars_ed = PSI.read_variables(results_ed)
    inj_pu = get_ptdf_bus_injections(ed_model)
    G, bus_lookup = build_graph_from_ptdf(inj_pu, ptdf; base_power = PSY.get_base_power(sys))
    add_hvdc_edges!(G, sys, res_vars_ed, bus_lookup; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches)
end

println("=== Scenario A: with circular flow (positive cost_sign, quadratic losses in ED) ===")
sys_a = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a, 1.0)
C_a = detect_cats_circular_flows_quad(sys_a)
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
C_b = detect_cats_circular_flows_quad(sys_b)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end
