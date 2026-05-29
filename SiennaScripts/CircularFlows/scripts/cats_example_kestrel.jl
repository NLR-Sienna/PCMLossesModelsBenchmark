# Circular-flow detection on the CATS system using a lossless PTDF UC model.
# Solved with HiGHS. Deploy at project root on Kestrel.
using Pkg
this_path = @__DIR__
Pkg.activate(this_path)
Pkg.instantiate()
kestrel_path = joinpath(this_path, "SiennaScripts", "CircularFlows")

using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerFlows
using PowerNetworkMatrices
using Graphs
using SimpleWeightedGraphs
using JuMP
using HiGHS
using Dates
using DataFrames
using Logging

include(joinpath(kestrel_path, "mapped_indices.jl"))
include(joinpath(kestrel_path, "circular_flows.jl"))
include(joinpath(this_path, "Systems", "CATS", "build_cats.jl"))
include(joinpath(this_path, "SiennaScripts", "build_models.jl"))
include(joinpath(this_path, "SiennaScripts", "utils.jl"))
include(joinpath(this_path, "SiennaScripts", "add_hvdc.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

cats_json = joinpath(this_path, "Systems", "CATS", "CATS_saved_reduced_sys.json")

highs_milp = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)

function get_ptdf_bus_injections(model)
    container = PSI.get_optimization_container(model)
    optimize!(PSI.get_jump_model(container))
    inj_expr = container.expressions[
        PSY.InfrastructureSystems.Optimization.ExpressionKey{
            PSI.ActivePowerBalance, PSY.ACBus}("")
    ]
    return JuMP.value.(inj_expr[:, 1]).data
end

function detect_cats_circular_flows(sys; time_step::Int = 1)
    ptdf = PTDF(sys)
    model = make_ptdf_model_without_losses(
        sys;
        device_models = CATS_UC_MODELS,
        optimizer     = highs_milp,
        ptdf          = ptdf,
        name          = "UC",
        ignore_pf     = true,
    )
    PSI.build!(model, output_dir = mktempdir())
    PSI.solve!(model)
    results = PSI.OptimizationProblemResults(model)
    res_vars = PSI.read_variables(results)
    inj_pu = get_ptdf_bus_injections(model)
    G, bus_lookup = build_graph_from_ptdf(inj_pu, ptdf; base_power = PSY.get_base_power(sys))
    add_hvdc_edges!(G, sys, res_vars, bus_lookup; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, bus_lookup, branches)
end

println("=== Scenario A: with circular flow (positive cost_sign, HVDC loop profitable) ===")
sys_a = build_cats_system(cats_json)
set_cats_renewable_costs!(sys_a, 1.0)
C_a = detect_cats_circular_flows(sys_a)
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
C_b = detect_cats_circular_flows(sys_b)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), min_flow=$(minimum(c.branch_flows)) MW")
end
