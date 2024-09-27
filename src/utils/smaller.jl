"""Sigtools.jl smaller()
This function is used to find the smaller of the two matrices in the n dimension.
This helps the same_length() function that makes matrices the same length.

Inputs:
    x::AbstractMatrix --> Data matrix 1
    y::AbstractMatrix --> Data matrix 2
    n::Int64 --> Dimension direction (range 1 or 2)

Outputs:
    x::AbstractMatrix

Sam Kramer

"""
function smaller(x::Matrix{T} where {T<:Real}, y::Matrix{T} where {T<:Real}, n::Int64)

  if size(x, n) < size(y, n)
    return x::AbstractMatrix
  elseif size(y, n) < size(x, n)
    return y::AbstractMatrix
  else
    return x::AbstractMatrix
  end

end

function smaller(x::Vector{T} where {T<:Real}, y::Vector{T} where {T<:Real})

  if length(x) < length(y)
    return x::Vector{T} where {T<:Real}
  elseif length(y) < length(x)
    return y::Vector{T} where {T<:Real}
  else
    return x::Vector{T} where {T<:Real}
  end

end
