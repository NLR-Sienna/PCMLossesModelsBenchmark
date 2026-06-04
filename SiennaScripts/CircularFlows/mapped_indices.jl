import Base: getindex, setindex!, size, isempty, axes, eachindex

"""
    MappedIndices <: AbstractArray{Int64,1}

A thin `AbstractArray` wrapper that lets you index a cycle vector of internal graph
node indices (as returned by `simplecycles`) and get back the corresponding PSY bus
numbers in one step (e.g. `bus_indices[cc]`).

Fields
------
- `indices`     : Values in internal-index space; the "contents" of the array.
                  Typically the sorted internal indices (1..N) produced by
                  `lookup_arrays`.
- `e2i_lookup`  : Dense lookup array; `e2i_lookup[external_bus_number]` returns the
                  corresponding internal index, or 0 if not present.
- `i2e_lookup`  : Dense lookup array; `i2e_lookup[internal_index]` returns the
                  corresponding PSY bus number (external ID).
"""
struct MappedIndices <: AbstractArray{Int64, 1}
    indices::Vector{Int64}
    e2i_lookup::Vector{Int64}
    i2e_lookup::Vector{Int64}
end

"""
    lookup_arrays(lookup::Dict{Int64,Int64}) -> (e2i_lookup, i2e_lookup, sorted_values)

Build the two dense lookup arrays used by `MappedIndices` from a
`bus_number => internal_index` Dict.

Returns
-------
- `e2i_lookup`     : `e2i_lookup[external_id]` gives the internal index, or 0 if
                     `external_id` is not a key in `lookup` (no `KeyError`).
- `i2e_lookup`     : `i2e_lookup[internal_index]` gives the external bus number.
- `sorted_values`  : Internal index values sorted in ascending order (used as the
                     `indices` field of `MappedIndices`).
"""
function lookup_arrays(lookup::Dict{Int64, Int64})
    k = collect(keys(lookup))
    v = collect(values(lookup))
    order = sortperm(v)
    if isempty(k) || isempty(v)
        e2i_lookup = Int64[]
        i2e_lookup = Int64[]
    else
        e2i_lookup = zeros(Int64, maximum(k))
        e2i_lookup[k[order]] .= v[order]
        i2e_lookup = zeros(Int64, maximum(v))
        i2e_lookup[v[order]] .= k[order]
    end
    return e2i_lookup, i2e_lookup, v[order]
end

"""
    MappedIndices(lookup::Dict{Int64,Int64})

Build a `MappedIndices` from a `bus_lookup` Dict as returned by
`build_graph_from_ptdf` or `build_graph_from_pf_aux_variables`
(mapping PSY bus number → internal graph index).
"""
function MappedIndices(lookup::Dict{Int64, Int64})
    e2i_lookup, i2e_lookup, v = lookup_arrays(lookup)
    return MappedIndices(v, e2i_lookup, i2e_lookup)
end

Base.size(A::MappedIndices) = size(A.indices)
Base.isempty(A::MappedIndices) = isempty(A.indices)
Base.axes(A::MappedIndices) = axes(A.indices)
Base.eachindex(A::MappedIndices) = eachindex(A.indices)

# `A[idx]` returns the external PSY bus number for the internal index stored at
# position `idx` in `A.indices`.
function Base.getindex(A::MappedIndices, idx)
    return A.i2e_lookup[A.indices[idx]]
end

function Base.setindex!(A::MappedIndices, value, idx)
    internal_index = A.e2i_lookup[value]
    internal_index == 0 && throw(ArgumentError("Value not found in lookup"))
    A.indices[idx] = internal_index
end

"""
    e2i(A::MappedIndices, value) -> Int64

Translate an external PSY bus number `value` to its internal graph index.
Throws `ArgumentError` if `value` is not present in the lookup.
"""
function e2i(A::MappedIndices, value)
    internal_index = A.e2i_lookup[value]
    internal_index == 0 && throw(ArgumentError("Value not found in lookup"))
    return internal_index
end

"""
    i2e(A::MappedIndices, idx) -> Int64

Translate a positional index `idx` into `A.indices` back to an external PSY bus number.
"""
function i2e(A::MappedIndices, idx)
    return A.i2e_lookup[A.indices[idx]]
end
