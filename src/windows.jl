#=  windows.jl

    This is the source code for the windows section of SigTools
    Equations taken from [1]Harris, "Use of Windows for Harmonic Analysis", 1976

    Sam Kramer, May 2024
=#

#NOTE: BASIC FUNCTIONS
"""
SigTools.jl windows.jl

This function will create a rectangular window function over length N where n_[0,N] = 0 and the rest are equal to 1.

Sam Kramer
"""
@inline function rectangular_window(N::Int64)

    duty_cycle = ones(Float64, 1, N - 2)
    window = append!(vec([0.0]), duty_cycle)
    window = append!(window, vec([0.0]))

    return window::Vector{Float64}
end

"""
SigTools.jl windows.jl

This function creates a triangular window function of length N from 0 to 1.

Sam Kramer
"""
@inline function triangle_window(N::Int64)

    duty_cycle = collect(LinRange(0.2, 1, floor(Int64, N / 2)))
    window = append!(duty_cycle, reverse(duty_cycle))

    if length(window) < N
        append!(window, 0.2)
    elseif length(window) > N
        pop!(window)
    end

    return window::Vector{Float64}
end

#NOTE: MATHEMATIC WINDOW FUNCTIONS
"""
SigTools.jl windows.jl

This function is will create a hann window around your selected data of length N. A Hann window is a sine windowing function. 

Inputs:
    N::Int64 == Length of vector output
    α::Int64 == Power of the sine function of the window

Sam Kramer
"""
@inline function hann(N::Int64, α::Int64)

    x = collect(LinRange(0, 1, N))
    window = sin.(π .* x) .^ α

    return window::Vector{Float64}
end

"""
SigTools.jl windows.jl

This function will create a Hamming window of length N.

Sam Kramer
"""
@inline function hamming(N::Int64)

    x = collect(LinRange(0, 1, N))
    window = 0.54 .- 0.46 .* cos.(2 .* π .* x)

    return window::Vector{Float64}
end

"""
SigTools.jl windows.jl

This function will create a Blackman window of length N.

Sam Kramer
"""
@inline function blackman(N::Int64)

    x = collect(LinRange(-0.5, 0.5, N))
    window = 0.42 .+ 0.50 .* cos.(2 .* π .* x) .+ 0.08 .* cos.(4 .* π .* x)

    return window::Vector{Float64}
end

#NOTE: CONSTRUCTED WINDOWS
"""
SigTools windows.jl

This function creates a constructed Riesz windowing function of length N.

Sam Kramer
"""
@inline function riesz(N::Int64)

    x = collect(LinRange(0, 0.5, floor(Int64, N / 2)))
    window = 1.0 .- abs.(2 .* x) .^ 2
    window = append!(reverse(window), window)

    if length(window) != N
        append!(window, window[end])
    end

    return window::Vector{Float64}
end

"""
SigTools.jl windows.jl

This function will generate a constructed Riemann function of length N.

Sam Kramer
"""
@inline function riemann(N::Int64)

    x = collect(LinRange(0, 0.5, floor(Int64, N / 2)))
    window = sin.(2 .* π .* x) ./ (2 .* π .* x)
    window = append!(reverse(window), window)

    if length(window) != N
        append!(window, window[end])
    end

    return window::Vector{Float64}
end

"""
SigTools.jl windows.jl

This function will generate a constructed Tukey winodw of length N with shape function of α.

Inputs:
    N::Int64 == Length of the output vector.
    α::Float64 == Shape vector of Tukey window function

Sam Kramer
"""
@inline function tukey(N::Int64, α::Float64)

    # --Check if α < N
    if α > 1.0
        error("Cannot have α be larger than 1. For an α = 0, the window returns a Rectangular window. For α = 1, the window returns a Hann window.")
    end

    # --Construct window
    x = collect(1:round(Int64, α * N / 2))
    window = sin.(π .* x ./ (α .* N)) .^ 2
    w = ones(Float64, 1, floor(Int64, N / 2 - α * N / 2))
    window = append!(window, w)
    second_half_window = reverse(window)
    window = append!(window, second_half_window)

    # --Check length
    if length(window) != N
        window = append!(window, zeros(Float64, 1, 10))
        window = window[1:N]
    end

    return window::Vector{Float64}
end

"""
SigTools.jl windows.jl

This is a function that will create a constructed Poisson windowing function of length N and shape α

Inputs:
    N::Int64 == Length of the output vector
    α::Float64 == Shape function of the window

Sam Kramer
"""
@inline function poisson(N::Int64, α::Float64)

    x = collect(1:round(N / 2))
    window = exp.(-α .* abs.(x) ./ (N / 2))
    reversed_window = reverse(window)
    window = append!(reversed_window, window)

    # --Check lengths
    if length(window) != N
        window = append!(window, zeros(Float64, 1, 10))
        window = window[1:N]
    end

    return window::Vector{Float64}
end
