# End-to-end validation of the CircularFlows module on the RTS test system.
# Two scenarios demonstrate correct circular flow detection:
#   Scenario A: renewable cost > 0  → HVDC loop is profitable → circular flow detected
#   Scenario B: renewable cost < 0  → renewables generate freely, no loop incentive → no circular flow

using Pkg
this_path = @__DIR__
#Pkg.activate(joinpath(this_path, "..", "..")) # Use this in your local machine
Pkg.activate(this_path) # Use this in Kestrel
Pkg.instantiate()
kestrel_path = joinpath(this_path, "SiennaScripts", "CircularFlows")


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
using Gurobi

#include(joinpath(this_path, "mapped_indices.jl"))
#include(joinpath(this_path, "circular_flows.jl"))
#include(joinpath(this_path, "../../Systems/RTS/build_rts.jl"))

include(joinpath(kestrel_path, "mapped_indices.jl"))
include(joinpath(kestrel_path, "circular_flows.jl"))
include(joinpath(kestrel_path, "../../Systems/RTS/build_rts.jl"))

const PSY = PowerSystems
const PSI = PowerSimulations

function detect_circular_flows(sys; time_step::Int=1)
    optimizer = PSI.optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
    #optimizer = PSI.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.01)
    model = PSI.DecisionModel(
        make_uc_template(), sys;
        optimizer=optimizer, name="UC", store_variable_names=true,
    )
    PSI.build!(model, output_dir=mktempdir())
    PSI.solve!(model)

    results = PSI.OptimizationProblemResults(model)
    res_vars = PSI.read_variables(results)
    data = PSI.get_power_flow_data(
        only(PSI.get_power_flow_evaluation_data(PSI.get_optimization_container(model)))
    )

    G = build_graph(data; time_step=time_step)
    add_hvdc_edges!(G, sys, res_vars, data; time_step=time_step)
    branches = collect(PSY.get_components(PSY.ACBranch, sys))
    return find_circular_flows(G, data, branches), model
end

# ---- Scenario A: negative renewable costs → circular flow expected ----
println("=== Scenario A: with circular flow (negative RE costs) ===")
sys_a = build_rts_system()
set_renewable_costs!(sys_a, 1.0) # Positive number is added negative to the objective function
C_a, model = detect_circular_flows(sys_a);
println("  Cycles found: $(length(C_a))")
for (i, c) in enumerate(C_a)
    println("  Cycle $i: bus_numbers=$(c.bus_numbers), branches=$(c.branches), min_flow=$(minimum(c.branch_flows)) MW")
end

# ---- Scenario B: HVDC disabled → no circular flow expected ----
# With no HVDC, there is no inter-area shortcut to create a directed loop in
# the AC network, so the circular-flow detector should return 0 cycles.
# Note: negative RE costs alone are insufficient to eliminate circular flows
# when HVDC is free (zero loss): the optimizer may still route power through
# HVDC+AC mesh loops. Disabling HVDC provides a clean negative control.
println("=== Scenario B: without circular flow (HVDC disabled) ===")
sys_b = build_rts_system()
set_renewable_costs!(sys_b, -1.0)
PSY.set_available!(only(PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys_b)), false)
C_b, model = detect_circular_flows(sys_b)
println("  Cycles found: $(length(C_b))")
for (i, c) in enumerate(C_b)
    println("  Cycle $i: bus_numbers=$(c.bus_numbers), branches=$(c.branches), min_flow=$(minimum(c.branch_flows)) MW")
end

