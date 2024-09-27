#=  convolutions.jl
    Basic defintions for convolutions for vectors and matrices.

    Sam Kramer
=#

using .SigTools
using FFTW
using Test

#-- 2D Convolutions
"""
"""
function convol(x::Matrix{T} where {T<:Real}, y::Matrix{T} where {T<:Real})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @inbounds conv2_xy = real.(ifft(f_x .* f_y))

    return conv2_xy::Matrix{Float64}
end

function convol(x::Matrix{Float64}, y::Matrix{Float64})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @inbounds conv2_xy = real.(ifft(f_x .* f_y))

    return conv2_xy::Matrix{Float64}
end

#-- 1D Convolutions
function convol(x::Vector{Float64}, y::Vector{Float64})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @inbounds conv_xy = real.(ifft(f_x .* f_y))

    return conv_xy::Vector{Float64}
end

function convol(x::Vector{T} where {T<:Real}, y::Vector{T} where {T<:Real})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @inbounds conv_xy = real.(ifft(f_x .* f_y))

    return conv_xy::Vector{T} where {T<:Real}
end

# --Deconvolutions
""" SigTools.jl deconvol
This is a deconvolution function that is defined as the deconvolution of y from x. This is done through the division of x by y in the frequency domain, which is defined by the convolution theorem. For image deconvolution, the noisy image should be x and y should be the noise model to deconvolve from the image.

Inputs:
    x::Vector{T} where {T<:Real} == Vector of data
    y::Vector{T} where {T<:Real} == Vector of data that will be deconvolved from x

Outputs:
    x::Vector{T} where {T<:Real} == Vector of data similar to x that has had y deconvolved from it

Sam Kramer
"""
function deconvol(x::Vector{T} where {T<:Real}, y::Vector{T} where {T<:Real})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @inbounds conv_xy = real.(ifft(f_x ./ f_y))

    return conv_xy::Vector{T} where {T<:Real}
end

function deconvol(x::Matrix{T} where {T<:Real}, y::Matrix{T} where {T<:Real})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @inbounds conv2_xy = real.(ifft(f_x ./ f_y))

    return conv2_xy::Matrix{Float64}
end

