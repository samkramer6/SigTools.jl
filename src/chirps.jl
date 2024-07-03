#=  chirps.jl
    This file is the description file for the chirps function for the SigTools package.

    Sam Kramer
    March 29th, 2024
=#

using .SigTools

###################################################################################################
#                                              Chirp                                              #
###################################################################################################

"""
SigTools.jl Chirp()
This is the function that describes the chrip function that is used to create chirp signals

Chirp(f_start, f_stop, T, cycles, type)
    Logarithmic --> Logarithmic frequency pattern
    Linear --> Linear frequency pattern
"""
function chirp(f_start, f_stop, T, cycles, type)

    # --Find fs needed
    fs = maximum([f_start, f_stop]) .* 1000

    # --Define time vector
    time = vec(collect(0:(1 / fs):(T / cycles - 1 / fs)))

    # --Define frequency vector
    if uppercase(type) == "LOGARITHMIC" || uppercase(type) == "LOG"

        f = vec(2 .* pi .* logrange(f_start, f_stop, length(time)))

    elseif uppercase(type) == "LINEAR"

        f = vec(2 .* pi .* LinRange(f_start, f_stop, length(time)))

    end

    # -- Define the chirp
    chirp_signal = cos.(f .* time)

    # --Define number of cycles
    for i in collect(1:cycles)
        append!(chirp_signal, chirp_signal)
    end

    # --Redefine time function
    time_out = vec(LinRange(0, T, length(chirp_signal)))

    return chirp_signal::Vector{Float64}, time_out
end

