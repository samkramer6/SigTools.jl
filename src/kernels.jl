#=  kernels.jl

    Function definitions for the different types of kernel functions that will be included within
    the julia package SigTools. Kernel functions are used for understanding data by transforming it
    either linearly or nonlinearly.

    Sam Kramer
    March 25th, 2024

=#

using .SigTools
using FFTW
using StatsBase
using LinearAlgebra
using Test

######################################################################################################
#                                           Linear Kernels                                           #
######################################################################################################

"""
    SigTools.jl linkernel()
    Linear kernel function

    Can compute the following:
        Absolute value
        Similarity Score
        Normalized Similarity Score

    Inputs:
        x::AbstractArray --> Data, can be matrix of size [N x M] or a vector
        y::AbstractArray --> Data, can be matrix of size [N x M] or a vector
        type::AbstractString --> String determining the type of kernel chosen

    Outputs:
        kx::AbstractArray or Float64 --> Output of Kernel Transformation, size depends on kernel selection


    Sam Kramer
    Methods:
"""
function linkernel(x::AbstractVecOrMat)

    kx = abs.(x)

    return kx::AbstractVecOrMat
end

function linkernel(x::AbstractVecOrMat, y::AbstractVecOrMat)

    (x, y) = same_length(x, y)
    @inbounds kx = abs.(x .+ y)

    return kx::AbstractVecOrMat
end

function linkernel(x::AbstractVecOrMat, y::AbstractVecOrMat, type::AbstractString)

    (x, y) = same_length(x, y)

    if uppercase(type) == "ABS" || uppercase(type) == "ABSOLUTE VALUE"

        kx = abs.(x .+ y)

    elseif uppercase(type) == "SIMSCORE" ||
           uppercase(type) == "SIMILARITY" ||
           uppercase(type) == "SIM"

        kx = dot(x, y)

    elseif uppercase(type) == "NORMSIMSCORE" || uppercase(type) == "NORM SIMILARITY"

        kx = dot(x, y) ./ length(x)

    else

        kx = x .+ y

    end

    return kx
end

######################################################################################################
#                                          Nonlinear Kernels                                         #
######################################################################################################


# --Default Kernels
"""
SigTools.jl kernel()
The default version of the nonlinear kernel function.

Computes following kernel functions with standard parameters:
    Gaussian: σ = 10
    Polynomial: n = 20
    Cosine Similarity: Euclidian Norm
    Laplacian: σ = 10

Inputs:
    x::AbstractArray --> Data, can be matrix of size [N x M] or a vector
    type::AbstractString --> String determining the type of kernel chosen

Outputs:
    kx::AbstractArray or Float64 --> Output of Kernel Transformation, size depends on kernel selection


Sam Kramer
Methods:
"""
function kernel(x::AbstractVecOrMat, type::AbstractString)      # 1 Input default NL-Kernel

    type = uppercase(type)

    if type == "GAUSSIAN" || type == "RBF"

        sigma = 10     # Width parameter
        @fastmath kx = exp.(-1 .* (x .^ 2) ./ (2 * sigma))

    elseif type == "POLY" || type == "POLYNOMIAL"

        n = 20
        @fastmath kx = (x .+ 1) .^ n

    elseif type == "COSSIM" || type == "SIMILARITY"

        x_norm = euclidian_norm(x)
        kx = sum(x) / x_norm

    elseif type == "LAPLACIAN"

        sigma = 10
        x_norm = euclidian_norm(x)
        kx = exp(-x_norm ./ 10)

    end

    return kx
end

function kernel(x::AbstractVecOrMat, y::AbstractVecOrMat, type::AbstractString)         # 2 Input default NL-Kernel

    (x, y) = same_length(x, y)
    type = uppercase(type)

    if type == "GAUSSIAN" || type == "RBF"

        sigma = 10
        @fastmath kx = exp.(-1 .* ((x .^ 2) + (y .^ 2)) ./ (2 * sigma))

    elseif type == "POLY" || type == "POLYNOMIAL"

        n = 20
        @fastmath kx = (x .* y .+ 1) .^ n

    elseif type == "COSSIM" || type == "SIMILARITY"

        x_norm = euclidian_norm(x::AbstractVecOrMat, y::AbstractVecOrMat)
        kx = sum(x) / x_norm

    elseif type == "LAPLACIAN"

        sigma = 10
        x_norm = euclidian_norm(x::AbstractVecOrMat, y::AbstractVecOrMat)
        kx = exp(-x_norm ./ 10)

    end

    return kx
end


# --Parameterized Kernel Functions

function kernel(
    x::Matrix{T} where {T<:Real},
    type::AbstractString,
    param::T where {T<:Real},
)

    type = uppercase(type)

    if type == "GAUSSIAN" || type == "RBF"

        sigma = param     # Width parameter
        @fastmath kx = exp.(-1 .* (x .^ 2) ./ (2 * sigma))

    elseif type == "POLY" || type == "POLYNOMIAL"

        n = param
        @fastmath kx = (x .+ 1) .^ n

    elseif type == "COSSIM" || type == "SIMILARITY"

        x_norm = euclidian_norm(x::AbstractVecOrMat)
        kx = sum(x) / x_norm

    elseif type == "LAPLACIAN"

        sigma = param
        x_norm = euclidian_norm(x)
        kx = exp(-1 .* x_norm ./ 10)

    end

    return kx
end

function kernel(x::AbstractVecOrMat, type::AbstractString, param::T where {T<:Real})

    type = uppercase(type)

    if type == "GAUSSIAN" || type == "RBF"

        sigma = param     # Width parameter
        @fastmath kx = exp.(-1 .* (x .^ 2) ./ (2 * sigma))

    elseif type == "POLY" || type == "POLYNOMIAL"

        n = param
        @fastmath kx = (x .+ 1) .^ n

    elseif type == "COSSIM" || type == "SIMILARITY"

        x_norm = euclidian_norm(x)
        kx = sum(x) / x_norm

    elseif type == "LAPLACIAN"

        sigma = param
        x_norm = euclidian_norm(x)
        kx = exp(-1 .* x_norm ./ 10)

    end

    return kx
end

function kernel(
    x::Matrix{T} where {T<:Real},
    y::Matrix{T} where {T<:Real},
    type::AbstractString,
    param::T where {T<:Real},
)

    (x, y) = same_length(x, y)
    type = uppercase(type)

    if type == "GAUSSIAN" || type == "RBF"

        sigma = param
        @fastmath kx = exp.(-1 .* (x .^ 2 - y .^ 2) ./ (2 * sigma))

    elseif type == "POLY" || type == "POLYNOMIAL"

        n = param
        @fastmath kx = (x .* y .+ 1) .^ n

    elseif type == "COSSIM" || type == "SIMILARITY"

        x_norm = euclidian_norm(x::AbstractVecOrMat, y::AbstractVecOrMat)
        kx = sum(x) / x_norm

    elseif type == "LAPLACIAN"

        sigma = param
        x_norm = euclidian_norm(x::AbstractArray, y::AbstractArray)
        kx = exp(-x_norm ./ 10)

    end

    return kx
end

function kernel(
    x::AbstractVecOrMat,
    y::AbstractVecOrMat,
    type::AbstractString,
    param::T where {T<:Real},
)

    (x, y) = same_length(x, y)
    type = uppercase(type)

    if type == "GAUSSIAN" || type == "RBF"

        sigma = param
        @fastmath kx = exp.(-1 .* (x .^ 2 - y .^ 2) ./ (2 * sigma))

    elseif type == "POLY" || type == "POLYNOMIAL"

        n = param
        @fastmath kx = (x .* y .+ 1) .^ n

    elseif type == "COSSIM" || type == "SIMILARITY"

        x_norm = euclidian_norm(x::AbstractArray, y::AbstractArray)
        kx = sum(x) / x_norm

    elseif type == "LAPLACIAN"

        sigma = param
        x_norm = euclidian_norm(x::AbstractArray, y::AbstractArray)
        kx = exp(-x_norm ./ 10)

    end

    return kx
end


######################################################################################################
#                                                Norms                                               #
######################################################################################################

"""
SigTools.jl euclidian_norm()
This is the function for calculating the Euclidian Norm of a data set

Inputs:
    x::AbstractArray --> Data, can be matrix of size [N x M] or a vector
    y::AbstractArray --> Data, can be matrix of size [N x M] or a vector


Sam Kramer
Methods:
"""
function euclidian_norm(x::AbstractVecOrMat)

    @inbounds norm = sqrt(sum(x .^ 2))

    return norm::Float64
end

function euclidian_norm(x::AbstractVecOrMat, y::AbstractVecOrMat)

    (x, y) = same_length(x, y)

    @inbounds norm = sqrt(sum((x .- y) .^ 2))

    return norm::Float64
end

"""
    SigTools.jl norm()
    This is the function to calculate a norm of any order

    Inputs:
        x::AbstractArray --> Data, can be matrix of size [N x M] or a vector
        y::AbstractArray --> Data, can be matrix of size [N x M] or a vector
        order::Int64 or Float64 --> Order of the norm, can be a Int or a Float to do partial norms


    Sam Kramer
    Methods:
"""
function norm(x::AbstractVecOrMat, order::T where {T<:Real})

    @fastmath norm_out = (sum(abs.(x) .^ order) .^ (1 / order))

    return norm_out::Float64
end

function norm(x::AbstractVecOrMat, y::AbstractVecOrMat, order::T where {T<:Real})

    (x, y) = same_length(x, y)

    @inbounds norm_out = (sum(abs.(x .- y) .^ order) .^ (1 / order))

    return norm_out::Float64
end

function norm(x::AbstractVecOrMat, y::AbstractVecOrMat, order::Float64)

    (x, y) = same_length(x, y)

    @inbounds norm_out = (sum(abs.(x .- y) .^ order) .^ (1 / order))

    return norm_out::T where {T<:Number}
end

"""
    SigTools.jl infnorm()
    This function is used to calculate the infinite norm of the data set.

    Inputs:
        x::AbstractArray --> Data, can be matrix of size [N x M] or a vector
        y::AbstractArray --> Data, can be matrix of size [N x M] or a vector


    Sam Kramer

    Methods:
"""
function infnorm(x::AbstractVecOrMat, y::AbstractVecOrMat)

    (x, y) = same_length(x, y)

    @inbounds norm_out = maximum(x .- y)

    return norm_out::T where {T<:Number}
end

function infnorm(x::AbstractVecOrMat)

    norm_out = maximum(x)

    return norm_out::T where {T<:Number}
end

