# Flow-cancelling transmission investment on the RTS GMLC system.
#
# Mirrors `script_5bus_flowcancelling.jl` but exercises the parallel-line
# (`-double_circuit`) and TapTransformer handling that RTS requires.
#
# Run:  julia --project=. script_rts_flowcancelling.jl
using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using Dates
using JuMP
using HiGHS
using Logging
using InfrastructureSystems
import PowerNetworkMatrices

include("Systems/RTS/build_rts.jl")
include("SiennaScripts/FlowCancelling/build_models.jl")
include("SiennaScripts/utils.jl")

const MILP = optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 1e-3)

"""
    run_rts_fc(; project_costs, choke_rating) -> (model, res, cand_names)

Build + solve the RTS flow-cancelling model. When `choke_rating` is given, the two series
segments that replace each candidate's original corridor have their rating set to that value
(p.u.), congesting the original path so the parallel candidate becomes useful — used by the
stressed scenario to force an investment.
"""
function run_rts_fc(; project_costs = Float64[], choke_rating = nothing)
    sys = build_rts_system()
    cand_names = add_candidate_lines_rts!(sys; n_candidates = 2, project_costs = project_costs)
    @info "RTS candidate lines added" cand_names

    if choke_rating !== nothing
        for cand in cand_names
            orig = replace(cand, "candidate_line_" => "")
            for seg in (orig * "_segment_1", orig * "_segment_2")
                line = get_component(Line, sys, seg)
                line === nothing || set_rating!(line, choke_rating)
            end
        end
    end

    model = build_model_with_flow_canceling_terms(
        sys;
        device_models = RTS_FC_DEVICE_MODELS,
        optimizer = MILP,
    )
    solve!(model)
    res = OptimizationProblemResults(model)
    return model, res, cand_names
end

"""Read z (build) and v (cancelling-flow) values from a solved model's container."""
function read_fc_results(model, cand_names)
    container = model.internal.container
    jump = container.JuMPmodel
    if !JuMP.has_values(jump)
        JuMP.optimize!(jump)
    end
    time_steps = PSI.get_time_steps(container)
    z = PSI.get_variable(container, BranchInvestmentVariable(), Line)
    v = PSI.get_variable(container, BranchCancellingFlowVariable(), Line)
    z_vals = Dict(c => JuMP.value(z[c]) for c in cand_names)
    v_max = Dict(c => maximum(t -> abs(JuMP.value(v[c, t])), time_steps) for c in cand_names)
    return z_vals, v_max
end

# ── Base scenario ──────────────────────────────────────────────────────────────
@info "=== RTS flow-cancelling: base scenario ==="
model, res, cand_names = run_rts_fc()

const VTOL = 1e-4

container = model.internal.container
status = JuMP.termination_status(container.JuMPmodel)
@info "termination status" status
@assert status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) "Model did not solve to optimality: $status"

z_vals, v_max = read_fc_results(model, cand_names)
println("\n=== Candidate investment results (base) ===")
for cand in cand_names
    println("  $cand : z = $(round(z_vals[cand], digits=4)),  max|v| = $(round(v_max[cand], digits=6))")
    if z_vals[cand] > 0.5
        @assert v_max[cand] <= VTOL "Built candidate $cand has nonzero v (max|v|=$(v_max[cand]))"
    end
end

# Sanity check: realized FC flows are finite for every monitored Line (incl. double_circuit).
if !JuMP.has_values(container.JuMPmodel)
    JuMP.optimize!(container.JuMPmodel)
end
fc_flow = PSI.get_expression(container, PTDFBranchFlowWithFC(), Line)
max_flow = 0.0
for name in axes(fc_flow, 1), t in PSI.get_time_steps(container)
    val = JuMP.value(fc_flow[name, t])
    @assert isfinite(val) "Non-finite FC flow on $name"
    global max_flow = max(max_flow, abs(val))
end
@info "Max |PTDFBranchFlowWithFC| over Lines (p.u.)" max_flow
@info "RTS flow-cancelling base scenario PASSED"

# ── Stressed scenario: free candidates + choked original corridors ──────────────
# Choking the segmented original path forces flow onto the parallel candidate, so at least
# one candidate should be built. Confirms the build/cancel mechanism is actually exercised.
@info "=== RTS flow-cancelling: stressed scenario (free candidates, choked corridors) ==="
model_s, _, cand_names_s = run_rts_fc(; project_costs = [0.0, 0.0], choke_rating = 0.05)
status_s = JuMP.termination_status(model_s.internal.container.JuMPmodel)
@assert status_s in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) "Stressed model not optimal: $status_s"

z_vals_s, v_max_s = read_fc_results(model_s, cand_names_s)
println("\n=== Candidate investment results (stressed) ===")
for cand in cand_names_s
    built = z_vals_s[cand] > 0.5
    println("  $cand : z = $(round(z_vals_s[cand], digits=4)),  max|v| = $(round(v_max_s[cand], digits=6)),  built=$built")
    if built
        @assert v_max_s[cand] <= VTOL "Built candidate $cand has nonzero v (max|v|=$(v_max_s[cand]))"
    end
end
@assert any(z_vals_s[c] > 0.5 for c in cand_names_s) "Stressed scenario built no candidate line"
@info "RTS flow-cancelling stressed scenario PASSED — at least one candidate built, v≈0 when built"

@info "ALL RTS FLOW-CANCELLING CHECKS PASSED"
