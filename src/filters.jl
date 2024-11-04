#   filters.jl

"""
SigTools.jl filters.jl

This is a filter data structure that will be immutable and used to describe the filter equation

Sam Kramer
"""
struct Filter
  order::Int64
  numerator::Vector{Float64}
  denominator::Vector{Float64}
end

"""
SigTools.jl filters.jl

This is a function that implements the filter data structure to compute a finite difference FIR filter of the data.

Inputs:
    filter_design::Filter ==
    data::Vector{Float64} ==

Outputs:
    filtered_data::Vector{Float64} == Filtered data

Sam Kramer
"""
function finite_diff_filt(filter_design::Filter, data::Vector{Float64})

  if length(filter_design.numerator) != filter_design.order + 1
    error("Order must be smaller than the number of coefficients. If Order 𝓞, is size N, then length of coefficients must be of size N + 1.")
  end
  if length(filter_design.numerator) != length(filter_design.denominator) + 1
    error("Denominator must be of length 1 smaller than the numerator")
  end

  A = filter_design.denominator
  B = filter_design.numerator
  # global filtered_data = vec(zeros(Float64, 1, length(data)))
  # global data = data
  exp = "$(B[1]) * data[i]"

  for n in collect(range(1, filter_design.order))
    exp *= " + $(B[1+n]) * data[i-$n] + $(A[n]) * filtered_data[i-$n]"
  end

  exp = Meta.parse(exp)

  # for i in collect(range(filter_design.order, length(data)))
  # filtered_data[i] = Meta.eval(exp)
  # end

  return exp
end

filter_design = Filter(2, vec([2 0.5 0.5]), vec([2 1]))
global data = vec(rand(Float64, 1, 10))
global output = finite_diff_filt(filter_design, data)
filtered_data = vec(zeros(Float64, 1, length(data)))
i = 2

@show output
Meta.eval(output)


