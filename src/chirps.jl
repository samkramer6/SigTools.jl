#TEST: chirp() function 

#=  chirps.jl
    This file is the description file for the chirps function for the SigTools package.

    Sam Kramer
    March 29th, 2024
=#

using .SigTools

###################################################################################################

struct Chirp
    signal::Vector{Float64}
    time::Vector{Float64}
end

###################################################################################################
#                                              Chirp                                              #
###################################################################################################

"""
SigTools.jl Chirp()
This is the function that describes the chrip function that is used to create chirp signals

chirp(f_start, f_stop, T, cycles, type)
    Logarithmic --> Logarithmic frequency pattern
    Linear --> Linear frequency pattern
"""
function chirp(
    f_start::T where {T<:Real},
    f_stop::T where {T<:Real},
    T::T where {T<:Real},
    type::String,
)

    # --Find fs needed
    const fs = maximum([f_start, f_stop]) .* 1000
    const cycles = 1;

    # --Define time vector
    const time = vec(collect(0:(1 / fs):(T / cycles - 1 / fs)))

    # --Define frequency vector
    if uppercase(type) == "LOGARITHMIC" || uppercase(type) == "LOG"
        f = vec(2 .* pi .* logrange(f_start, f_stop, length(time)))
    elseif uppercase(type) == "LINEAR"
        f = vec(2 .* pi .* LinRange(f_start, f_stop, length(time)))
    end

    # -- Define the chirp and time vector
    chirp_signal::Vector{Float64} = cos.(f .* time)
    time_out::Vector{Float64} = vec(LinRange(0, T, length(chirp_signal)))

    return Chirp(chirp_signal, time_out)
end

function chirp(
    f_start::T where {T<:Real},
    f_stop::T where {T<:Real},
    T::T where {T<:Real},
    cycles::Int64,
    type::String,
)

    # --Find fs needed
    const fs = maximum([f_start, f_stop]) .* 1000

    # --Define time vector
    const time = vec(collect(0:(1 / fs):(T / cycles - 1 / fs)))

    # --Define frequency vector
    if uppercase(type) == "LOGARITHMIC" || uppercase(type) == "LOG"
        f = vec(2 .* pi .* logrange(f_start, f_stop, length(time)))
    elseif uppercase(type) == "LINEAR"
        f = vec(2 .* pi .* LinRange(f_start, f_stop, length(time)))
    end

    # -- Define the chirp and time vector
    chirp_signal::Vector{Float64} = cos.(f .* time)
        Iterators.flatten(Iterators.repeated(append!(chirp_signal, chirp_signal), cycles));
    time_out::Vector{Float64} = vec(LinRange(0, T, length(chirp_signal)))

    return Chirp(chirp_signal, time_out)
end

function chirp(
    f_start::T where {T<:Real},
    f_stop::T where {T<:Real},
    T::T where {T<:Real},
    cycles::Int64,
    type::String,
    fs::Int64,
)
    # --Define time vector
    const time = vec(collect(0:(1 / fs):(T / cycles - 1 / fs)))

    # --Define frequency vector
    if uppercase(type) == "LOGARITHMIC" || uppercase(type) == "LOG"
        f = vec(2 .* pi .* logrange(f_start, f_stop, length(time)))
    elseif uppercase(type) == "LINEAR"
        f = vec(2 .* pi .* LinRange(f_start, f_stop, length(time)))
    end

    # -- Define the chirp and time vector
    chirp_signal::Vector{Float64} = cos.(f .* time)
        Iterators.flatten(Iterators.repeated(append!(chirp_signal, chirp_signal), cycles));
    time_out::Vector{Float64} = vec(LinRange(0, T, length(chirp_signal)))

    return Chirp(chirp_signal, time_out)
end

