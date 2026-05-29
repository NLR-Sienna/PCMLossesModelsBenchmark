import SimpleWeightedGraphs as SWG
import Graphs as GR
using PowerSystems
using DataFrames
using Dates

const PSY = PowerSystems

struct CircularFlow
    buses::Vector{Int64}        # internal graph node indices
    bus_numbers::Vector{Int64}  # PSY bus numbers
    branches::Vector{String}    # branch names in cycle
    branch_flows::Vector{Float64}  # flow magnitudes in MW
end

# ── Private API helpers ────────────────────────────────────────────────────────

# Returns Dict{bus_number => internal_index} from PSY5 PowerFlowData.
function _get_bus_lookup(data)::Dict{Int64, Int64}
    arcs = data.power_network_matrix.arc_admittance_from_to.axes[1]
    all_bus_numbers = sort!(unique(vcat(first.(arcs), last.(arcs))))
    return Dict{Int64, Int64}(bn => i for (i, bn) in enumerate(all_bus_numbers))
end

# Returns arc iterator as (from_internal_idx, to_internal_idx) pairs.
function _get_arc_iter(data, bus_lookup::Dict{Int64, Int64})
    arcs = data.power_network_matrix.arc_admittance_from_to.axes[1]
    return Tuple{Int64, Int64}[(bus_lookup[first(a)], bus_lookup[last(a)]) for a in arcs]
end

# Returns (flows_from_to, flows_to_from) in MW.
function _get_flows(data, time_step::Int, flow_type::Symbol, base_power::Float64)
    if flow_type == :active
        return (data.arc_active_power_flow_from_to[:, time_step] .* base_power,
                data.arc_active_power_flow_to_from[:, time_step] .* base_power)
    else
        return (data.arc_reactive_power_flow_from_to[:, time_step] .* base_power,
                data.arc_reactive_power_flow_to_from[:, time_step] .* base_power)
    end
end

# Shared HVDC edge insertion: adds one directed edge per HVDC line.
# get_flow is a closure that returns the from-to flow in per-unit for a given HVDC component.
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
internal sequential indices (1..N_buses); use `_get_bus_lookup` to map back to
bus numbers. Works with both the old and new PowerFlowData API variants.
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
to real AC branch loops. Uses `simplecycles` from Graphs.jl.
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

Use this when the AC power flow evaluation fails to converge; the PTDF flows
are exact for the lossless DC model the optimizer solved against.
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

Overload for graphs built with `build_graph_from_ptdf`. Accepts `bus_lookup`
(Dict{bus_number => internal_index}) directly instead of PowerFlowData.
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
