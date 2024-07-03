#=  logrange.jl

    This is the function that will be used to define the logrange function in SigTools.jl

    Sam Kramer
    March 24th, 2024
=#

###################################################################################################
#                                              logrange                                           #
###################################################################################################

"""
  SigTools.jl   logrange()
  This is the function definition for the logrange() function that defines a logarithmic range from
  input 1 to input 2 with n steps.

  Inputs:
    x1::AbstractInt
    x2::AbstractFloat -->
    n::Int64 --> Number of steps between x1 and x2

  Outputs:
    log_x::Vec{AbstractFloat} --> Vector of the logarithmic range from x1 to x2.

  Sam Kramer

  Methods:
"""
function logrange(x1::Float64, x2::Float64, n::Int)

    # --Compute the power needed
    stop = log10(x2)
    start = log10(x1)

    # --Compute the range needed
    @fastmath log_x = 10 .^ (range(start, stop, length = n))

    return log_x::Vector{Float64}
end

function logrange(x1::Int, x2::Int, n::Int)

    # --Compute power needed
    stop = log10(x2)
    start = log10(x1)

    # --Computer the range needed
    @fastmath log_x = 10 .^ (range(start, stop, length = n))

    return log_x::Vector{Float64}
end
