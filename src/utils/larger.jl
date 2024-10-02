"""Sigtools.jl larger()
This function is used to find the larger of the two matrices in the n dimension.
This helps the same_length() function that makes matrices the same length.

Inputs:
    x::AbstractMatrix --> Data matrix 1
    y::AbstractMatrix --> Data matrix 2
    n::Int64 --> Dimension direction (range 1 or 2)

Outputs:
    x::AbstractMatrix

Sam Kramer
"""
function larger(x::Vector{T} where {T<:Real}, y::Vector{T} where {T<:Real})

    if length(x) > length(y)
        return x::AbstractVector
    elseif length(y) > length(x)
        return y::AbstractVector
    else
        return x::AbstractVector
    end

end

function larger(x::Matrix{T} where {T<:Real}, y::Matrix{T} where {T<:Real}, n::Int64)

    if size(x, n) > size(y, n)
        return x::AbstractArray
    elseif size(y, n) > size(x, n)
        return y::AbstractArray
    else
        return x::AbstractArray
    end

end

