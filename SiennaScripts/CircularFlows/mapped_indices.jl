import Base: getindex, setindex!, size, isempty, axes, eachindex

struct MappedIndices <: AbstractArray{Int64, 1}
    indices::Vector{Int64}
    e2i_lookup::Vector{Int64}
    i2e_lookup::Vector{Int64}
end

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

function MappedIndices(lookup::Dict{Int64, Int64})
    e2i_lookup, i2e_lookup, v = lookup_arrays(lookup)
    return MappedIndices(v, e2i_lookup, i2e_lookup)
end

Base.size(A::MappedIndices) = size(A.indices)
Base.isempty(A::MappedIndices) = isempty(A.indices)
Base.axes(A::MappedIndices) = axes(A.indices)
Base.eachindex(A::MappedIndices) = eachindex(A.indices)

function Base.getindex(A::MappedIndices, idx)
    return A.i2e_lookup[A.indices[idx]]
end

function Base.setindex!(A::MappedIndices, value, idx)
    internal_index = A.e2i_lookup[value]
    internal_index == 0 && throw(ArgumentError("Value not found in lookup"))
    A.indices[idx] = internal_index
end

function e2i(A::MappedIndices, value)
    internal_index = A.e2i_lookup[value]
    internal_index == 0 && throw(ArgumentError("Value not found in lookup"))
    return internal_index
end

function i2e(A::MappedIndices, idx)
    return A.i2e_lookup[A.indices[idx]]
end
