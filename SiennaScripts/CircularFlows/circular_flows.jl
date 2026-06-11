import SimpleWeightedGraphs as SWG
import Graphs as GR
using PowerSystems
using DataFrames
using Dates

const PSY = PowerSystems

# Both naming conventions that PSI uses across formulation versions.
const _HVDC_FLOW_KEYS = (
    "FlowActivePowerVariable__TwoTerminalGenericHVDCLine",
    "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine",
)

"""
    read_hvdc_flow_variables(res [, read_fn]) -> Dict{String, DataFrame}

Return a `Dict` containing whichever HVDC flow variable key(s) are present in `res`.
Missing keys are silently skipped, so the Dict may be empty (which `add_hvdc_edges!`
handles gracefully with an @info notice).

`read_fn` is the PSI function used to read a single variable by name:
- `read_variable`          for `OptimizationProblemResults`  (default)
- `read_realized_variable` for `SimulationProblemResults`

Avoids calling `read_variables(res)` for *all* stored variables, which would fail
when post-build variables such as `LineLossTotalApproximation__System` are present in
the container but were never allocated in the serialised results.
"""
function read_hvdc_flow_variables(res, read_fn::Function = read_variable)::Dict{String, DataFrame}
    out = Dict{String, DataFrame}()
    for k in _HVDC_FLOW_KEYS
        try
            out[k] = read_fn(res, k)
        catch
        end
    end
    return out
end

"""
    CircularFlow

A single detected loop-flow (circular flow) in the power network.

Fields
------
- `buses`        : Internal graph node indices (1..N) for each bus in the cycle.
                   These are NOT PSY bus numbers; use `bus_numbers` for human-readable IDs.
- `bus_numbers`  : PSY bus numbers corresponding to each node in `buses`.
- `branches`     : Names of the AC branches (and any HVDC links) that form the cycle,
                   in the same traversal order as `buses`.
- `branch_flows` : MW flow magnitude on each branch in `branches` (always positive;
                   direction is encoded by the edge orientation in the graph).
"""
struct CircularFlow
    buses::Vector{Int64}        # internal graph node indices
    bus_numbers::Vector{Int64}  # PSY bus numbers
    branches::Vector{String}    # branch names in cycle
    branch_flows::Vector{Float64}  # flow magnitudes in MW
end

# ── Private API helpers ────────────────────────────────────────────────────────

# Returns Dict{bus_number => internal_index} built from the union of all from/to bus numbers
# appearing in the arc admittance axes of a PowerFlowData object.
function _get_bus_lookup(data)::Dict{Int64, Int64}
    arcs = data.power_network_matrix.arc_admittance_from_to.axes[1]
    all_bus_numbers = sort!(unique(vcat(first.(arcs), last.(arcs))))
    return Dict{Int64, Int64}(bn => i for (i, bn) in enumerate(all_bus_numbers))
end

# Returns a Vector of (from_internal_idx, to_internal_idx) tuples, one per arc,
# with bus numbers translated to internal graph indices via bus_lookup.
function _get_arc_iter(data, bus_lookup::Dict{Int64, Int64})
    arcs = data.power_network_matrix.arc_admittance_from_to.axes[1]
    return Tuple{Int64, Int64}[(bus_lookup[first(a)], bus_lookup[last(a)]) for a in arcs]
end

# Returns (flows_from_to, flows_to_from) in MW for the given time_step and flow_type
# (:active or :reactive), scaling per-unit PowerFlowData values by base_power.
function _get_flows(data, time_step::Int, flow_type::Symbol, base_power::Float64)
    if flow_type == :active
        return (data.arc_active_power_flow_from_to[:, time_step] .* base_power,
                data.arc_active_power_flow_to_from[:, time_step] .* base_power)
    else
        return (data.arc_reactive_power_flow_from_to[:, time_step] .* base_power,
                data.arc_reactive_power_flow_to_from[:, time_step] .* base_power)
    end
end

# Shared HVDC edge insertion: iterates over all TwoTerminalGenericHVDCLine components in sys
# and adds one directed edge per line. `get_flow` must be a closure accepting an HVDC component
# and returning its from-to flow in per-unit (positive = power flows from→to bus).
function _add_hvdc_edges_inner!(
    G::SWG.SimpleWeightedDiGraph,
    sys::PSY.System,
    bus_lookup::Dict{Int64, Int64},
    get_flow::Function,
)
    for h in PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys)
        f = bus_lookup[PSY.get_number(PSY.get_from_bus(h))]
        t = bus_lookup[PSY.get_number(PSY.get_to_bus(h))]
        flow_ft = get_flow(h)
        if flow_ft > 0.0
            GR.add_edge!(G, f, t, flow_ft)
        elseif flow_ft < 0.0
            GR.add_edge!(G, t, f, -flow_ft)
        # zero-flow: skip (no net HVDC power → no edge in graph)
        end
    end
end

# ── Public API ─────────────────────────────────────────────────────────────────

"""
    build_graph(data; base_power, time_step, flow_type) -> SWG.SimpleWeightedDiGraph

Build a directed flow graph from AC power flow data.

Each arc is oriented in the direction of net positive flow. Graph nodes are
**internal sequential indices (1..N_buses)** — they are NOT PSY bus numbers.
Callers must use `_get_bus_lookup(data)` to obtain the mapping from PSY bus
numbers to these internal indices (and vice versa) before interpreting node IDs.
Edge weights are in MW (i.e., `base_power × per-unit flow`).
Works with both the old and new PowerFlowData API variants.
"""
function build_graph(
    data;
    base_power::Float64 = 100.0,
    time_step::Int = 1,
    flow_type::Symbol = :active,
)
    flow_type in (:active, :reactive) ||
        throw(ArgumentError("flow_type must be :active or :reactive"))

    bus_lookup = _get_bus_lookup(data)
    arc_iter = _get_arc_iter(data, bus_lookup)
    flows_ft, flows_tf = _get_flows(data, time_step, flow_type, base_power)

    src = Vector{Int64}()
    dst = Vector{Int64}()
    w = Vector{Float64}()

    for (i, (f, t)) in enumerate(arc_iter)
        flow_ft = flows_ft[i]
        flow_tf = flows_tf[i]
        direction_ft = true
        if sign(flow_ft) == sign(flow_tf)
            abs(flow_ft) < abs(flow_tf) && (direction_ft = false)
        elseif sign(flow_ft) != 1
            direction_ft = false
        end
        if direction_ft
            push!(src, f); push!(dst, t); push!(w, flow_ft)
        else
            push!(src, t); push!(dst, f); push!(w, flow_tf)
        end
    end

    return SWG.SimpleWeightedDiGraph(src, dst, w)
end

"""
    add_hvdc_edges!(G, sys, res_vars::Dict, data; time_step::Int)

RTS-style: add HVDC edges from a PSI `read_variables` result Dict.
Handles both `FlowActivePowerVariable__...` and `FlowActivePowerFromToVariable__...` keys.
`res_vars` is typically produced by `read_hvdc_flow_variables(res)`, which
selects only the HVDC-relevant keys from the PSI results container.
"""
function add_hvdc_edges!(
    G::SWG.SimpleWeightedDiGraph,
    sys::PSY.System,
    res_vars::Dict,
    data;
    time_step::Int = 1,
)
    bus_lookup = _get_bus_lookup(data)
    if haskey(res_vars, "FlowActivePowerVariable__TwoTerminalGenericHVDCLine")
        hvdc_df = res_vars["FlowActivePowerVariable__TwoTerminalGenericHVDCLine"]
    elseif haskey(res_vars, "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine")
        hvdc_df = res_vars["FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine"]
    else
        @info "No HVDC flow variable found in results (HVDC may be disabled); skipping HVDC edges"
        return
    end
    hvdc_df isa Pair && (hvdc_df = last(hvdc_df))  # unwrap Pair{DateTime, DataFrame} if needed
    # Detect DataFrame format: PSI can return either
    #   (a) long format:  columns = ["DateTime", "name", "value"]  (newer PSI versions)
    #   (b) wide format:  columns = ["DateTime", component_name1, component_name2, ...]
    is_long_format = "name" in names(hvdc_df) && "value" in names(hvdc_df)
    if is_long_format
        # Long format: filter by component name, then select the time_step-th occurrence
        function _get_hvdc_flow_long(h)
            comp_name = PSY.get_name(h)
            rows = hvdc_df[hvdc_df.name .== comp_name, :value]
            return Float64(rows[time_step])
        end
        _add_hvdc_edges_inner!(G, sys, bus_lookup, _get_hvdc_flow_long)
    else
        # Wide format: each component has its own column
        # Identify numeric columns to exclude the DateTime column
        hvdc_cols = [n for n in names(hvdc_df) if eltype(hvdc_df[!, n]) <: Number]
        hvdc_names_in_sys = [PSY.get_name(h) for h in PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys)]
        # Exact name match first; fall back to positional order
        function _get_hvdc_col_wide(h_name)
            if h_name in hvdc_cols
                return h_name
            else
                idx = findfirst(==(h_name), hvdc_names_in_sys)
                return isnothing(idx) ? hvdc_cols[1] : hvdc_cols[idx]
            end
        end
        _add_hvdc_edges_inner!(G, sys, bus_lookup, h -> Float64(hvdc_df[time_step, _get_hvdc_col_wide(PSY.get_name(h))]))
    end
end

"""
    add_hvdc_edges!(G, sys, hvdc_flows_ft::DataFrame, data; time_step::DateTime)

CATS-style: add HVDC edges from an explicit long-format DataFrame with columns
`DateTime`, `name`, `value`. Used when the HVDC result is extracted separately
from the simulation results (as in the CATS workflow).
`time_step` is a `DateTime` (not an Int) because it is used to filter rows of
the long-format DataFrame by their `DateTime` column value.
"""
function add_hvdc_edges!(
    G::SWG.SimpleWeightedDiGraph,
    sys::PSY.System,
    hvdc_flows_ft::DataFrame,
    data;
    time_step::DateTime,
)
    bus_lookup = _get_bus_lookup(data)
    _add_hvdc_edges_inner!(G, sys, bus_lookup, h -> begin
        vals = hvdc_flows_ft[
            (hvdc_flows_ft.DateTime .== time_step) .& (hvdc_flows_ft.name .== PSY.get_name(h)),
            :value,
        ]
        vals[1]
    end)
end

"""
    find_circular_flows(G, data, branches) -> Vector{CircularFlow}

Enumerate all simple directed cycles in G and return those that correspond
to real AC branch loops. Uses `simplecycles` (Johnson's algorithm,
O((V+E)(C+1)) where C is the number of simple cycles) from Graphs.jl.
For large, highly meshed grids the number of cycles can be enormous;
callers may want to pre-filter the graph or limit system size before calling.
"""
function find_circular_flows(
    G::SWG.SimpleWeightedDiGraph,
    data,
    branches::Vector{<:PSY.ACBranch},
)
    bus_lookup = _get_bus_lookup(data)
    bus_indices = MappedIndices(bus_lookup)
    result = CircularFlow[]
    for cc in GR.simplecycles(G)
        c_branches, branch_flows = _find_connected_branches(G, bus_lookup, branches, cc)
        push!(result, CircularFlow(cc, bus_indices[cc], PSY.get_name.(c_branches), branch_flows))
    end
    return result
end

"""
    build_graph_from_ptdf(bus_injection_pu, ptdf; base_power) -> (G, bus_lookup)

Build a directed flow graph from DC branch flows computed via PTDF.
`bus_injection_pu` must be a plain vector of per-unit bus injections in the same
bus ordering as ptdf (i.e., position k corresponds to ptdf[arc, k]).

Use this when AC power flow was not run (e.g., `ignore_pf_ed = true`), or as a
cross-check against the AC-flow-based graph. Branch flows are computed as
`PTDF × bus_injection_pu × base_power`, which is exact for the lossless DC model
the optimizer solved against.
"""
function build_graph_from_ptdf(
    bus_injection_pu::AbstractVector{Float64},
    ptdf;
    base_power::Float64 = 100.0,
)
    arc_ax = ptdf.axes[2]   # (from_bus, to_bus) arc tuples, length = n_arcs
    num_arcs = length(arc_ax)
    num_buses = length(bus_injection_pu)

    all_bus_numbers = sort!(unique(vcat(first.(arc_ax), last.(arc_ax))))
    bus_lookup = Dict{Int64, Int64}(bn => i for (i, bn) in enumerate(all_bus_numbers))

    src = Vector{Int64}()
    dst = Vector{Int64}()
    w   = Vector{Float64}()

    for (j, (from_bus, to_bus)) in enumerate(arc_ax)
        flow_mw = sum(ptdf[j, k] * bus_injection_pu[k] for k in 1:num_buses) * base_power
        f = bus_lookup[from_bus]
        t = bus_lookup[to_bus]
        if flow_mw > 0
            push!(src, f); push!(dst, t); push!(w, flow_mw)
        elseif flow_mw < 0
            push!(src, t); push!(dst, f); push!(w, -flow_mw)
        end
    end

    return SWG.SimpleWeightedDiGraph(src, dst, w), bus_lookup
end

"""
    build_graph_from_pf_aux_variables(res_ed, sys; time_step) -> (G, bus_lookup)

Build a directed flow graph from actual AC power flow branch flows stored in PSI
aux variables after an ED solve with `power_flow_evaluation` enabled.
These aux variables are only present when the model was built with
`ignore_pf_ed = false`; calling this function on results from an
`ignore_pf_ed = true` solve will error.

Reads `PowerFlowBranchActivePowerFromTo/ToFrom__Line`,
`PowerFlowBranchActivePowerFromTo/ToFrom__TapTransformer` (RTS), and
`PowerFlowBranchActivePowerFromTo/ToFrom__Transformer2W` (CATS) — all in MW,
already scaled — and applies the same direction-selection logic as `build_graph`.
Both transformer reads are wrapped in try/catch so systems that lack a given
transformer type degrade gracefully (empty DataFrames are substituted).
Branches whose names are not found in `sys` (e.g., filtered or renamed components)
are silently skipped.

Returns `(G, bus_lookup)` with the same signature as `build_graph_from_ptdf`
so `add_hvdc_edges!(G, sys, res_vars, bus_lookup)` and
`find_circular_flows(G, bus_lookup, branches)` can be called unchanged.
"""
function build_graph_from_pf_aux_variables(
    res_ed,
    sys::PSY.System;
    time_step::Int = 1,
)
    line_ft_df = read_realized_aux_variable(
        res_ed, "PowerFlowBranchActivePowerFromTo__Line"; table_format = TableFormat.WIDE)
    line_tf_df = read_realized_aux_variable(
        res_ed, "PowerFlowBranchActivePowerToFrom__Line"; table_format = TableFormat.WIDE)
    local tap_ft_df, tap_tf_df
    try
        tap_ft_df = read_realized_aux_variable(
            res_ed, "PowerFlowBranchActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)
        tap_tf_df = read_realized_aux_variable(
            res_ed, "PowerFlowBranchActivePowerToFrom__TapTransformer"; table_format = TableFormat.WIDE)
    catch
        tap_ft_df = DataFrame()
        tap_tf_df = DataFrame()
    end
    local t2w_ft_df, t2w_tf_df
    try
        t2w_ft_df = read_realized_aux_variable(
            res_ed, "PowerFlowBranchActivePowerFromTo__Transformer2W"; table_format = TableFormat.WIDE)
        t2w_tf_df = read_realized_aux_variable(
            res_ed, "PowerFlowBranchActivePowerToFrom__Transformer2W"; table_format = TableFormat.WIDE)
    catch
        t2w_ft_df = DataFrame()
        t2w_tf_df = DataFrame()
    end

    line_names = names(line_ft_df)[2:end]
    tap_names  = isempty(tap_ft_df) ? String[] : names(tap_ft_df)[2:end]
    t2w_names  = isempty(t2w_ft_df) ? String[] : names(t2w_ft_df)[2:end]

    all_bus_numbers = sort!(unique(Int64[
        [PSY.get_number(PSY.get_from_bus(PSY.get_component(PSY.Line, sys, n)))           for n in line_names];
        [PSY.get_number(PSY.get_to_bus(PSY.get_component(PSY.Line, sys, n)))             for n in line_names];
        [PSY.get_number(PSY.get_from_bus(PSY.get_component(PSY.TapTransformer, sys, n))) for n in tap_names];
        [PSY.get_number(PSY.get_to_bus(PSY.get_component(PSY.TapTransformer, sys, n)))   for n in tap_names];
        [PSY.get_number(PSY.get_from_bus(PSY.get_component(PSY.Transformer2W, sys, n)))  for n in t2w_names];
        [PSY.get_number(PSY.get_to_bus(PSY.get_component(PSY.Transformer2W, sys, n)))    for n in t2w_names];
    ]))
    bus_lookup = Dict{Int64, Int64}(bn => i for (i, bn) in enumerate(all_bus_numbers))

    src = Vector{Int64}()
    dst = Vector{Int64}()
    w   = Vector{Float64}()

    function _push_edge!(flow_ft, flow_tf, from_no, to_no)
        (flow_ft == 0.0 && flow_tf == 0.0) && return
        f = bus_lookup[from_no]
        t = bus_lookup[to_no]
        direction_ft = true
        if sign(flow_ft) == sign(flow_tf)
            abs(flow_ft) < abs(flow_tf) && (direction_ft = false)
        elseif sign(flow_ft) != 1
            direction_ft = false
        end
        if direction_ft
            push!(src, f); push!(dst, t); push!(w, flow_ft)
        else
            push!(src, t); push!(dst, f); push!(w, flow_tf)
        end
    end

    for name in line_names
        branch = PSY.get_component(PSY.Line, sys, name)
        isnothing(branch) && continue
        _push_edge!(
            Float64(line_ft_df[time_step, name]),
            Float64(line_tf_df[time_step, name]),
            PSY.get_number(PSY.get_from_bus(branch)),
            PSY.get_number(PSY.get_to_bus(branch)),
        )
    end

    for name in tap_names
        branch = PSY.get_component(PSY.TapTransformer, sys, name)
        isnothing(branch) && continue
        _push_edge!(
            Float64(tap_ft_df[time_step, name]),
            Float64(tap_tf_df[time_step, name]),
            PSY.get_number(PSY.get_from_bus(branch)),
            PSY.get_number(PSY.get_to_bus(branch)),
        )
    end

    for name in t2w_names
        branch = PSY.get_component(PSY.Transformer2W, sys, name)
        isnothing(branch) && continue
        _push_edge!(
            Float64(t2w_ft_df[time_step, name]),
            Float64(t2w_tf_df[time_step, name]),
            PSY.get_number(PSY.get_from_bus(branch)),
            PSY.get_number(PSY.get_to_bus(branch)),
        )
    end

    return SWG.SimpleWeightedDiGraph(src, dst, w), bus_lookup
end

"""
    build_graph_from_acopf_variables(res_ed, sys; time_step) -> (G, bus_lookup)

Build a directed flow graph from `ACPPowerModel` optimization variables.
Use this when the ED was built with `ignore_pf_ed = true` (no post-solve power flow).
For the post-solve-PF mode use `build_graph_from_pf_aux_variables` instead.

Reads `FlowActivePowerFromToVariable__Line` and `FlowActivePowerToFromVariable__Line`
(per-unit), plus optional transformer variables for CATS (`__Transformer2W`) and
RTS (`__TapTransformer`), and scales all flows to MW using `sys` base power.

Returns `(G, bus_lookup)` with the same signature as `build_graph_from_ptdf`
so `add_hvdc_edges!` and `find_circular_flows` can be called unchanged.
"""
function build_graph_from_acopf_variables(
    res_ed,
    sys::PSY.System;
    time_step::Int = 1,
)
    base_power = PSY.get_base_power(sys)
    ptdf = PTDF(sys)  # for bus_lookup construction only; flows come from variables, not PTDF
    PNM.populate_branch_maps_by_type!(ptdf.network_reduction_data)

    line_ft_df = read_realized_variable(
        res_ed, "FlowActivePowerFromToVariable__Line"; table_format = TableFormat.WIDE)
    line_tf_df = read_realized_variable(
        res_ed, "FlowActivePowerToFromVariable__Line"; table_format = TableFormat.WIDE)

    local t2w_ft_df, t2w_tf_df
    try
        t2w_ft_df = read_realized_variable(
            res_ed, "FlowActivePowerFromToVariable__Transformer2W"; table_format = TableFormat.WIDE)
        t2w_tf_df = read_realized_variable(
            res_ed, "FlowActivePowerToFromVariable__Transformer2W"; table_format = TableFormat.WIDE)
    catch
        t2w_ft_df = DataFrame()
        t2w_tf_df = DataFrame()
    end

    local tap_ft_df, tap_tf_df
    try
        tap_ft_df = read_realized_variable(
            res_ed, "FlowActivePowerFromToVariable__TapTransformer"; table_format = TableFormat.WIDE)
        tap_tf_df = read_realized_variable(
            res_ed, "FlowActivePowerToFromVariable__TapTransformer"; table_format = TableFormat.WIDE)
    catch
        tap_ft_df = DataFrame()
        tap_tf_df = DataFrame()
    end

    line_names = names(line_ft_df)[2:end]
    t2w_names  = isempty(t2w_ft_df) ? String[] : names(t2w_ft_df)[2:end]
    tap_names  = isempty(tap_ft_df) ? String[] : names(tap_ft_df)[2:end]

    # Build reverse map: aggregated variable name → [individual line names].
    # component_to_reduction_name_map[Line] maps individual → aggregated;
    # reverse it so aggregated variable names can be resolved to bus numbers.
    nrd = ptdf.network_reduction_data
    agg_to_lines = Dict{String, Vector{String}}()
    if haskey(nrd.component_to_reduction_name_map, Line)
        for (ind_name, agg_name) in nrd.component_to_reduction_name_map[Line]
            push!(get!(agg_to_lines, agg_name, String[]), ind_name)
        end
    end

    # Resolve a Line variable name (direct or aggregated) to (from_bus_no, to_bus_no).
    # Returns nothing if unresolvable. Parallel-circuit variable names (e.g.
    # "A33-double_circuit") are not PSY components; look up a constituent line instead.
    function _line_endpoints(name)
        branch = PSY.get_component(PSY.Line, sys, name)
        !isnothing(branch) && return (PSY.get_number(PSY.get_from_bus(branch)),
                                       PSY.get_number(PSY.get_to_bus(branch)))
        if haskey(agg_to_lines, name)
            rep = PSY.get_component(PSY.Line, sys, first(agg_to_lines[name]))
            !isnothing(rep) && return (PSY.get_number(PSY.get_from_bus(rep)),
                                        PSY.get_number(PSY.get_to_bus(rep)))
        end
        return nothing
    end

    line_endpoints   = filter(!isnothing, [_line_endpoints(n) for n in line_names])
    t2w_branches     = filter(!isnothing, [PSY.get_component(PSY.Transformer2W, sys, n) for n in t2w_names])
    t2w_from_num_bus = [PSY.get_number(PSY.get_from_bus(b)) for b in t2w_branches]
    t2w_to_num_bus   = [PSY.get_number(PSY.get_to_bus(b))   for b in t2w_branches]
    tap_branches     = filter(!isnothing, [PSY.get_component(PSY.TapTransformer, sys, n) for n in tap_names])
    tap_from_num_bus = [PSY.get_number(PSY.get_from_bus(b)) for b in tap_branches]
    tap_to_num_bus   = [PSY.get_number(PSY.get_to_bus(b))   for b in tap_branches]
    all_bus_numbers = sort!(unique(Int64[
        first.(line_endpoints);
        last.(line_endpoints);
        t2w_from_num_bus;
        t2w_to_num_bus;
        tap_from_num_bus;
        tap_to_num_bus;
    ]))
    bus_lookup = Dict{Int64, Int64}(bn => i for (i, bn) in enumerate(all_bus_numbers))

    src = Vector{Int64}()
    dst = Vector{Int64}()
    w   = Vector{Float64}()

    function _push_edge!(flow_ft_pu, flow_tf_pu, from_no, to_no)
        flow_ft = Float64(flow_ft_pu) * base_power
        flow_tf = Float64(flow_tf_pu) * base_power
        (flow_ft == 0.0 && flow_tf == 0.0) && return
        f = bus_lookup[from_no]
        t = bus_lookup[to_no]
        direction_ft = true
        if sign(flow_ft) == sign(flow_tf)
            abs(flow_ft) < abs(flow_tf) && (direction_ft = false)
        elseif sign(flow_ft) != 1
            direction_ft = false
        end
        if direction_ft
            push!(src, f); push!(dst, t); push!(w, flow_ft)
        else
            push!(src, t); push!(dst, f); push!(w, flow_tf)
        end
    end

    for name in line_names
        endpoints = _line_endpoints(name)
        isnothing(endpoints) && continue
        from_no, to_no = endpoints
        _push_edge!(
            line_ft_df[time_step, name],
            line_tf_df[time_step, name],
            from_no,
            to_no,
        )
    end

    for name in t2w_names
        branch = PSY.get_component(PSY.Transformer2W, sys, name)
        isnothing(branch) && continue
        _push_edge!(
            t2w_ft_df[time_step, name],
            t2w_tf_df[time_step, name],
            PSY.get_number(PSY.get_from_bus(branch)),
            PSY.get_number(PSY.get_to_bus(branch)),
        )
    end

    for name in tap_names
        branch = PSY.get_component(PSY.TapTransformer, sys, name)
        isnothing(branch) && continue
        _push_edge!(
            tap_ft_df[time_step, name],
            tap_tf_df[time_step, name],
            PSY.get_number(PSY.get_from_bus(branch)),
            PSY.get_number(PSY.get_to_bus(branch)),
        )
    end

    return SWG.SimpleWeightedDiGraph(src, dst, w), bus_lookup
end

"""
    add_hvdc_edges!(G, sys, res_vars::Dict, bus_lookup::Dict; time_step)

Overload for graphs built with `build_graph_from_ptdf`. Accepts `bus_lookup`
(Dict{bus_number => internal_index}) directly instead of PowerFlowData.
"""
function add_hvdc_edges!(
    G::SWG.SimpleWeightedDiGraph,
    sys::PSY.System,
    res_vars::Dict,
    bus_lookup::Dict{Int64, Int64};
    time_step::Int = 1,
)
    if haskey(res_vars, "FlowActivePowerVariable__TwoTerminalGenericHVDCLine")
        hvdc_df = res_vars["FlowActivePowerVariable__TwoTerminalGenericHVDCLine"]
    elseif haskey(res_vars, "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine")
        hvdc_df = res_vars["FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine"]
    else
        @info "No HVDC flow variable found in results (HVDC may be disabled); skipping HVDC edges"
        return
    end
    hvdc_df isa Pair && (hvdc_df = last(hvdc_df))
    is_long_format = "name" in names(hvdc_df) && "value" in names(hvdc_df)
    if is_long_format
        _add_hvdc_edges_inner!(G, sys, bus_lookup,
            h -> Float64(hvdc_df[hvdc_df.name .== PSY.get_name(h), :value][time_step]))
    else
        hvdc_cols = [n for n in names(hvdc_df) if eltype(hvdc_df[!, n]) <: Number]
        hvdc_names_in_sys = [PSY.get_name(h) for h in PSY.get_components(PSY.TwoTerminalGenericHVDCLine, sys)]
        function _col(h_name)
            h_name in hvdc_cols ? h_name : let
                idx = findfirst(==(h_name), hvdc_names_in_sys)
                isnothing(idx) ? hvdc_cols[1] : hvdc_cols[idx]
            end
        end
        _add_hvdc_edges_inner!(G, sys, bus_lookup,
            h -> Float64(hvdc_df[time_step, _col(PSY.get_name(h))]))
    end
end

"""
    find_circular_flows(G, bus_lookup, branches) -> Vector{CircularFlow}

Overload for graphs built with `build_graph_from_ptdf` or
`build_graph_from_pf_aux_variables`. Accepts `bus_lookup`
(Dict{bus_number => internal_index}) directly instead of PowerFlowData.
Uses `simplecycles` (Johnson's algorithm, O((V+E)(C+1)) where C is the number
of simple cycles); callers should guard against running this on large meshed grids
where C may be very large.
"""
function find_circular_flows(
    G::SWG.SimpleWeightedDiGraph,
    bus_lookup::Dict{Int64, Int64},
    branches::Vector{<:PSY.ACBranch},
)
    bus_indices = MappedIndices(bus_lookup)
    result = CircularFlow[]
    for cc in GR.simplecycles(G)
        c_branches, branch_flows = _find_connected_branches(G, bus_lookup, branches, cc)
        push!(result, CircularFlow(cc, bus_indices[cc], PSY.get_name.(c_branches), branch_flows))
    end
    return result
end

# Private: find AC branches whose endpoints are both in the given cycle.
function _find_connected_branches(
    G::SWG.SimpleWeightedDiGraph,
    bus_lookup::Dict{Int64, Int64},
    branches::Vector{<:PSY.ACBranch},
    cycle::Vector{Int64},
)
    connected = PSY.ACBranch[]
    flows = Float64[]
    n = length(cycle)
    cycle_edges = Set{Tuple{Int64, Int64}}()
    for i in 1:n
        push!(cycle_edges, (cycle[i], cycle[mod1(i + 1, n)]))
        push!(cycle_edges, (cycle[mod1(i + 1, n)], cycle[i]))
    end
    for b in branches
        f = bus_lookup[PSY.get_number(PSY.get_from_bus(b))]
        t = bus_lookup[PSY.get_number(PSY.get_to_bus(b))]
        if (f, t) in cycle_edges
            push!(connected, b)
            w = G.weights[f, t]
            w == 0.0 && (w = G.weights[t, f])
            push!(flows, w)
        end
    end
    return connected, flows
end
