import Base: getindex, setindex!, size, isempty, axes, eachindex

# Define the struct
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

function MappedIndices(indices::Vector{Int64}, lookup::Dict{Int64, Int64})
    e2i_lookup, i2e_lookup, v = lookup_arrays(lookup)
    return MappedIndices(sort(indices), e2i_lookup, i2e_lookup)
end

# Implement the AbstractArray interface
Base.size(A::MappedIndices) = size(A.indices)
Base.isempty(A::MappedIndices) = isempty(A.indices)
Base.axes(A::MappedIndices) = axes(A.indices)
Base.eachindex(A::MappedIndices) = eachindex(A.indices)

# Overload getindex
function getindex(A::MappedIndices, idx)
    return A.i2e_lookup[A.indices[idx]]
end

# Overload setindex!
function setindex!(A::MappedIndices, value, idx)
    internal_index = A.e2i_lookup[value]
    if any(internal_index .== 0)
        throw(ArgumentError("Value not found in lookup"))
    end
    A.indices[idx] = internal_index
end

function e2i(A::MappedIndices, value)
    internal_index = A.e2i_lookup[value]
    if any(internal_index .== 0)
        throw(ArgumentError("Value not found in lookup"))
    end
    return internal_index
end

function i2e(A::MappedIndices, idx)
    return A.i2e_lookup[A.indices[idx]]
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Example usage
    indices = [3, 1, 4]
    lookup = Dict(11 => 1, 12 => 2, 13 => 3, 14 => 4, 15 => 5)
    c = MappedIndices(lookup)

    c1 = MappedIndices(indices, lookup)

    # Accessing the array using user-facing values
    println(c[1])  # Output: 11
    println(c.indices[1])  # Output: 1

    println(c1[1])  # Output: 11
    println(c1.indices[1])  # Output: 1

    # Setting a value using user-facing values and translating the value
    c[1] = 12
    println(c[1])  # Output: 13
    println(c.indices[1])  # Output: 3

    c1[1] = 13
    println(c1[1])  # Output: 13
    println(c1.indices[1])  # Output: 3
end

