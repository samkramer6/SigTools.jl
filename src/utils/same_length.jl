""" SigTools.jl   same_length()
This is a utility function that can be used to create two vectors or matrices which are the same length. It zero pads the two vectors so that they are the same size. Two matrices A and B are [m x n] and [p x r] where m > p and r > n, the two outputs will both be of size [m x r].

Inputs:
  x -- Input Vector/Matrix
  y -- Input Vector/Matrix 

Outputs:
  Tuple{x, y} -- Tuple output of x and y where they are both the same size

Sam Kramer
"""
function same_length(x::AbstractVector, y::AbstractVector)

  if length(x) == length(y)
    return x, y
  end

  leng_diff = abs.(length(x) - length(y))
  pad = vec(zeros(typeof(x[1]), 1, leng_diff))

  if length(x) > length(y)

    @inbounds append!(y, pad)

  elseif length(y) > length(x)

    @inbounds append!(x, pad)

  end

  return x::AbstractVector, y::AbstractVector
end

function same_length(x::Vector{Float64}, y::Vector{Float64})

  if length(x) == length(y)
    return x, y
  end

  leng_diff = abs.(length(x) - length(y))
  pad = vec(zeros(Float64, 1, leng_diff))

  if length(x) > length(y)

    @inbounds append!(y, pad)

  elseif length(y) > length(x)

    @inbounds append!(x, pad)

  end

  return x::Vector{Float64}, y::Vector{Float64}
end

########################################################################
#                            2D same_length()                          #
########################################################################

function same_length(x::Matrix{T} where {T<:Real}, y::Matrix{T} where {T<:Real})

  # --Test
  if size(x) == size(y)
    return x::AbstractMatrix, y::AbstractMatrix
  end

  # --Plan Height
  height_plan = size(larger(x, y, 1), 1)
  width_plan = size(larger(x, y, 2), 2)

  x_new = zeros(height_plan, width_plan)
  y_new = zeros(height_plan, width_plan)

  x_new[1:size(x, 1), 1:size(x, 2)] = x
  y_new[1:size(y, 1), 1:size(y, 2)] = y

  return x_new::Matrix{Float64}, y_new::Matrix{Float64}
end

function same_length(x::Matrix{Float64}, y::Matrix{Float64})

  # --Test
  if size(x) == size(y)
    return x::Matrix{Float64}, y::Matrix{Float64}
  end

  # --Plan Height
  height_plan = size(larger(x, y, 1), 1)
  width_plan = size(larger(x, y, 2), 2)

  # --Heigth Padding
  height_diff = abs.(height_plan - size(smaller(x, y, 1), 1))
  height_pad = zeros(height_diff, size(smaller(x, y, 1), 2))

  if length(height_pad) != 0

    @inbounds a = vcat(smaller(x, y, 1), height_pad)
    b = larger(x, y, 1)

  elseif length(height_pad) == 0

    a = x
    b = y

  end

  # --Width Padding
  width_diff = abs.(width_plan - size(smaller(x, y, 2), 2))
  width_pad = zeros(width_diff, size(smaller(x, y, 2), 1))

  if length(width_pad) != 0

    pad = zeros(height_plan, width_diff)
    @inbounds c = hcat(smaller(a, b, 2), pad)
    d = larger(a, b, 2)


  elseif length(width_pad) == 0

    c = a
    d = b

  end

  if size(c) == size(d)
    return c::Matrix{Float64}, d::Matrix{Float64}
  else
    (c, d) == same_length(c, d)
  end
end

function same_length(x::Matrix{Int64}, y::Matrix{Int64})

  x = convert(Matrix{Float64}, x)
  y = convert(Matrix{Float64}, y)

  (x, y) = same_length(x, y)

  return x::Matrix{Float64}, y::Matrix{Float64}
end
