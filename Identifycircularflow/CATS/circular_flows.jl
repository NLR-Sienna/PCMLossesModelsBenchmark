using SimpleWeightedGraphs
using Graphs
using SparseArrays
using PowerSystems
using DataFrames
using PowerFlows
import PowerFlows: PowerFlowData

const PSY = PowerSystems

# using Sandbox


struct CircularFlow
    buses::Vector{Int64}
    bus_numbers::Vector{Int64}
    branches::Vector{String}
    branch_flows::Vector{Float64}
end


function build_graph(data; base_power::Float64=100., time_step::Int64=1, flow_type::Symbol=:active)

    src = Vector{Int64}()
    dst = Vector{Int64}()
    w = Vector{Float64}()

    if flow_type == :active
        flows_ft = data.arc_active_power_flow_from_to[:, time_step] .* base_power
        flows_tf = data.arc_active_power_flow_to_from[:, time_step] .* base_power
    elseif flow_type == :reactive
        flows_ft = data.arc_reactive_power_flow_from_to[:, time_step] .* base_power
        flows_tf = data.arc_reactive_power_flow_to_from[:, time_step] .* base_power
    else
        throw(ArgumentError("flow_type must be either :active or :reactive"))
    end

    for (i, (f, t)) in enumerate(zip(first.(data.power_network_matrix.arc_admittance_from_to.axes[1]), last.(data.power_network_matrix.arc_admittance_from_to.axes[1])))
        flow_ft = flows_ft[i]
        flow_tf = flows_tf[i]
        direction_ft = true
        if sign(flow_ft) == sign(flow_tf)
            if abs(flow_ft) < abs(flow_tf)
                direction_ft = false
            end
        elseif sign(flow_ft) != 1
            direction_ft = false
        end

        if direction_ft
            push!(src, f)
            push!(dst, t)
            push!(w, flow_ft)
        else
            push!(src, t)
            push!(dst, f)
            push!(w, flow_tf)
        end
    end

    G = SimpleWeightedDiGraph(src, dst, w)
    return G
end

#function add_hvdc_edges!(G::SimpleWeightedDiGraph, sys::System, res_vars::Dict, data::PowerFlowData; time_step::DateTime)
function add_hvdc_edges!(G::SimpleWeightedDiGraph, sys::System, hvdc_flows_ft::DataFrame, data::PowerFlowData; time_step::DateTime)
#=
    if haskey(res_vars, "FlowActivePowerVariable__TwoTerminalGenericHVDCLine")  # when using ACPPowerModel
        # withdrawal at from bus, injection at to bus
        hvdc_flows_ft = res_vars["FlowActivePowerVariable__TwoTerminalGenericHVDCLine"]
    elseif haskey(res_vars, "FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine")  # when using PTDF or DCPPowerModel
        # withdrawal at from bus, injection at to bus
        hvdc_flows_ft = res_vars["FlowActivePowerFromToVariable__TwoTerminalGenericHVDCLine"]
    else
        # no HVDC flows available
        return
    end
=#

    if typeof(hvdc_flows_ft) != DataFrame
        _, hvdc_flows_ft = first(hvdc_flows_ft)
    end
    hvdc = collect(get_components(PSY.TwoTerminalGenericHVDCLine, sys))
    for h in hvdc
        f = get_number(PSY.get_from_bus(h))
        t = get_number(PSY.get_to_bus(h))
        flow_ft = hvdc_flows_ft[(hvdc_flows_ft.DateTime .== time_step) .&(hvdc_flows_ft.name .== get_name(h)), :value]
        flow_tf = -hvdc_flows_ft[(hvdc_flows_ft.DateTime .== time_step) .&(hvdc_flows_ft.name .== get_name(h)), :value]
        direction_ft = true
        if sign(flow_ft[1]) == sign(flow_tf[1])
            if abs(flow_ft[1]) < abs(flow_tf[1])
                direction_ft = false
            end
        elseif sign(flow_ft[1]) != 1
            direction_ft = false
        end

        if direction_ft
            add_edge!(G, f, t, flow_ft[1])
            # @show "Adding HVDC edge from $f to $t with flow $flow_ft"
        else
            add_edge!(G, t, f, flow_tf[1])
            # @show "Adding HVDC edge from $t to $f with flow $flow_tf"
        end
    end
end



function find_circular_flows(G::SimpleWeightedDiGraph, data::PowerFlowData, branches::Vector{PSY.ACBranch})
    first_vals = data.power_network_matrix.arc_admittance_from_to.axes[2]
    second_vals = data.power_network_matrix.arc_admittance_from_to.axes[2]
    need_data = Dict(k => v for (k, v) in zip(first_vals, second_vals))
    bus_indices = MappedIndices(need_data)
    cycles = simplecycles(G)
    cycle_objects = Vector{CircularFlow}()
    for cc in cycles
        c_branches, branch_flows = find_connected_branches(G, data, branches, bus_indices, cc)
        C = CircularFlow(cc, bus_indices[cc], get_name.(c_branches), branch_flows)
        push!(cycle_objects, C)
    end
    return cycle_objects
end

function find_connected_branches(G::SimpleWeightedDiGraph, data::PowerFlowData, branches::Vector{PSY.ACBranch}, bus_indices::MappedIndices, cycle::Vector{Int64})
    #bus_lookup = data.bus_lookup
    connected_branches = Vector{PSY.ACBranch}()
    branch_flows = Vector{Float64}()
    for b in branches
        f, t = get_number(PSY.get_from_bus(b)), get_number(PSY.get_to_bus(b))
        if (f in cycle) && (t in cycle)
            push!(connected_branches, b)
            w = G.weights[f, t]
            if w == 0.0
                w = G.weights[t, f]
            end
            push!(branch_flows, w)
        end
    end
    return connected_branches, branch_flows
end
