using FFTW
using Distributions
using DSP

"""
SigTools.jl additive white gaussian noise awgn() function

This function is meant to simulate additive white Gaussian noise in data sets. It is mean-zeroed and
of specified length N.

Inputs:
    N::Int64 -- Number of data points needed to be returned by the data set
    A::Float64 -- Absolute value of the amplitude of the maximum values

Sam Kramer
"""
@inline function awgn(N::Int64, A::T where {T<:Real})

    mean = 0.0
    std = 1

    noise = rand(Normal(mean, std), 1, N)
    noise = sqrt(2) .* A .* vec(noise) ./ maximum(noise)

    return noise::Vector{Float64}
end

"""
SigTools.jl additive white gaussian noise awgn() function

This function is meant to simulate additive white Gaussian noise in data sets. It is not mean-zeroed and
centered at σ of specified length N.

Inputs:
    σ::Real -- The center of the noise.
    N::Int64 -- Number of data points needed to be returned by the data set
    A::Float64 -- Absolute value of the amplitude of the maximum values


Sam Kramer
"""
function white_noise(σ::T where {T<:Real}, N::Int64, A::T where {T<:Real})

    noise = awgn(N, A) .+ σ

    return noise::Vector{Float64}
end


"""
SigTools.jl Colored Noise processes

These group of functions will be meant to emulate colored noise processes

These follow a white noise + filter design

This follows an autoregressive process (ARP) filter which is of order 63

Sam Kramer
"""
function red_noise(N::Int64, A::T where {T<:Real})

    # --Create white noise data vector
    @inline noise = awgn(N, A)

    # --Filter that data to make colored noise
    a = vec(zeros(1, 63))
    a[1] = 1
    i = collect(range(2, length(a)))
    @inbounds a[i] .= (0.3 .+ i .- 1) .* (a[i .- 1] ./ i)
    filter_design = DSP.PolynomialRatio(vec(a), vec([1]))

    # --Color the noise
    red_noise = DSP.filt(filter_design, noise)

    return red_noise::Vector{Float64}
end

""" SigTools.jl noise.jl

This is a function that generates pink noise using a white noise + filter process. The filter is a 3 pole and 3 zero filter used to color the noise.

Input:
    N::Int64 == length of the vector output
    A::T where {T<:Real} == Amplitude constant

Output:
    pink_noise::Vector{Float64}

Sam Kramer
"""
function deep_red_noise(N::Int64, A::T where {T<:Real})

    # --Create white noise data vector
    noise = awgn(N, A)

    # --Design Pinkening filter
    zero = [0.98443604, 0.83392334, 0.07568359]
    pole = [0.99572754, 0.94790649, 0.53567505]
    filter_design = DSP.ZeroPoleGain(zero, pole, 1.0)

    # --filter data
    pink_noise = DSP.filt(filter_design, noise)

    return pink_noise::Vector{Float64}
end

""" SigTools.jl noise.jl

This is a function that generates pink noise using a Voss generator algorithm. The filter is a 3 pole and 3 zero filter used to color the noise.

Input:
    N::Int64 == length of the vector output
    A::T where {T<:Real} == Amplitude constant

Output:
    pink_noise::Vector{Float64}

Sam Kramer
"""
function voss_noise(N::Int64, A::T where {T<:Real})

    @inline row_0 = awgn(N, 1.0)
    row_1 = vec(zeros(Float64, 1, N))
    row_2 = vec(zeros(Float64, 1, N))
    row_3 = vec(zeros(Float64, 1, N))
    row_4 = vec(zeros(Float64, 1, N))
    rand_start::Int64 = rand(1:(N / 2))

    for i in collect(range(1, N - 17))

        if rem(i, 2) == 0
            @inbounds row_1[i:(i + 1)] = row_0[rand_start:(rand_start + 1)]
        elseif rem(i, 4) == 0
            @inbounds row_2[i:(i + 3)] = row_0[rand_start:(rand_start + 3)]
        elseif rem(i, 8) == 0
            @inbounds row_3[i:(i + 7)] = row_0[rand_start:(rand_start + 7)]
        elseif rem(i, 16) == 0
            @inbounds row_4[i:(i + 15)] = row_0[rand_start:(rand_start + 15)]
        end

    end

    noise_output = row_0 .+ row_1 .+ row_2 .+ row_3 .+ row_4
    noise_output = A .* sqrt.(mean.(noise_output .^ 2 ./ maximum(noise_output)))
    noise_output = noise_output .- mean(noise_output)

    return noise_output::Vector{Float64}
end


