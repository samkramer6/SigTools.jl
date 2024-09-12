#=  noise.jl

    This is going to be the noise script that is meant to help with the simulation of noise types
    for DSP purposes.

    Sam Kramer
=#

# --Using statements
using FFTW
using Distributions
using DSP

###################################################################################################
#                                             White Gaussian News                                 #
###################################################################################################

"""
SigTools.jl additive white gaussian noise awgn() function

This function is meant to simulate additive white Gaussian noise in data sets. It is mean-zeroed and
of specified length N.

Inputs:
    N::Int64 -- Number of data points needed to be returned by the data set
    A::Float64 -- Absolute value of the amplitude of the maximum values

"""
function awgn(N::Int64, A::T where {T<:Real})

    mean = 0.0
    std = 1

    noise = rand(Normal(mean, std), 1, N)
    noise = sqrt(2) .* A .* vec(noise) ./ maximum(noise)

    return noise::Vector{Float64}
end

function white_noise(σ::T where {T<:Real}, N::Int64, A::T where {T<:Real})

    noise = awgn(N, A) .+ σ

    return noise::Vector{Float64}
end

###################################################################################################
#                                        Colored Noise                                            #
###################################################################################################

"""
SigTools.jl Colored Noise processes

These group of functions will be meant to emulate colored noise processes

These follow a white noise + filter design

This follows an autoregressive process (ARP) filter which is of order 63
       _______ 
      |       |
----> | noise | ----> filter ----> gain
      |_______|

"""
function red_noise(N::Int64, A::T where {T<:Real})

    # --Create white noise data vector
    noise = awgn(N, A)

    # --Filter that data to make colored noise
    a = vec(zeros(1, 255))
    a[1] = 1
    i = collect(range(2, length(a)))
    a[i] .= (0.3 .+ i .- 1) .* (a[i.-1] ./ i)
    filter_design = DSP.PolynomialRatio(vec(a), vec([1]))

    # --Color the noise
    colored_noise = DSP.filt(filter_design, noise)

    return colored_noise::Vector{Float64}
end

###################################################################################################
#
###################################################################################################



