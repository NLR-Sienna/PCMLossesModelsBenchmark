# Run from the scripts/ directory on a local machine:
#   julia SiennaScripts/CircularFlows/scripts/rts_example_local.jl
using Pkg
this_path = @__DIR__  # = .../SiennaScripts/CircularFlows/scripts
Pkg.activate(joinpath(this_path, "..", "..", ".."))
Pkg.instantiate()

using PowerSystemCaseBuilder
using PowerSystems
using PowerSimulations
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using Graphs
using SimpleWeightedGraphs
using HiGHS
using Dates
using DataFrames

include(joinpath(this_path, "..", "mapped_indices.jl"))
include(joinpath(this_path, "..", "circular_flows.jl"))
include(joinpath(this_path, "..", "..", "..", "Systems", "RTS", "build_rts.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

function detect_circular_flows(sys; time_step::Int = 1)
    optimizer = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
    model = PSI.DecisionModel(
        make_uc_template(), sys;
        optimizer = optimizer, name = "UC", store_variable_names = true,
    )
    PSI.build!(model, output_dir = mktempdir())
    PSI.solve!(model)
    results = PSI.OptimizationProblemResults(model)
    res_vars = PSI.read_variables(results)
    data = PSI.get_power_flow_data(
        only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(model)))
    )
    G = build_graph(data; time_step = time_step)
    add_hvdc_edges!(G, sys, res_vars, data; time_step = time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, data, branches), model
end

println("=== Scenario A: with circular flow (positive RE cost_sign) ===")
sys_a = build_rts_system()
set_renewable_costs!(sys_a, 1.0)
C_a, _ = detect_circular_flows(sys_a)
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: buses=$(c.bus_numbers), branches=$(c.branches), min_flow=$(minimum(c.branch_flows)) MW")
end

println("=== Scenario B: without circular flow (HVDC disabled) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0)
PSY.set_available!(only(PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys_b)), false)
C_b, _ = detect_circular_flows(sys_b)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: buses=$(c.bus_numbers), branches=$(c.branches), min_flow=$(minimum(c.branch_flows)) MW")
end
