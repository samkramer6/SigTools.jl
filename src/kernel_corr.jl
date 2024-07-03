#=  kernel_corr.jl

    This file is going to be the defintion of the kernel correlation functions

    Sam Kramer
    March 26th, 2024
=#

using .SigTools
using FFTW
using Test

################################################################################################### 
#                                             KXCORR                                              #
###################################################################################################

"""
    SigTools.jl kxcorr()
    Kernel auto correlation and kernel cross-correlation function defintions

    inputs:
        x::AbstractArray -->
        y::AbstractArray -->
        kernel_type::AbstractString -->
        param::Float64 --> 

    outputs:
        corr_out::AbstractArray --> 

    Sam Kramer
    Methods:
"""
function kxcorr(x::Vector{T} where T <: Real, kernel_type::AbstractString, param::T where T <: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end  

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        corr_out = corr(kx);

    return corr_out::Vector{T} where T <: Real
end

function kxcorr(x::Vector{Float64}, kernel_type::AbstractString, param::T where T <: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end 

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        corr_out = corr(kx);

    return corr_out::Vector{Float64}
end

function kxcorr(x::Vector{Float64}, y::Vector{Float64}, kernel_type::AbstractString, param::T where T <: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end 

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        @inline ky = kernel(y, kernel_type, param);
        corr_out = corr(kx, ky);

    return corr_out::Vector{Float64}
end

function kxcorr(x::Vector{T} where T <: Real, y::Vector{T} where T <: Real, kernel_type::AbstractString, param::T where T <: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end 

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        @inline ky = kernel(y, kernel_type, param);
        corr_out = corr(kx, ky);

    return corr_out::Vector{T} where T <: Real
end

###################################################################################################
#                                   2D KXCORR for Matrices                                        #
###################################################################################################

function kxcorr(x::AbstractMatrix, kernel_type::AbstractString, param::T where T <: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end  

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        corr_out = corr(kx); 

    return corr_out::AbstractMatrix
end

function kxcorr(x::Matrix{Float64}, kernel_type::AbstractString, param::T where T<: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end 

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        corr_out = corr(kx);

    return corr_out::Matrix{T} where T <: Real
end

function kxcorr(x::Matrix{T} where T <: Real, y::Matrix{T} where T <: Real, kernel_type::AbstractString, param::T where T <: Real)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end

    # --Same length()
        (x, y) == same_length(x, y);

    # --Find Corr
        @inline kx = kernel(x, kernel_type, param);
        @inline ky = kernel(y, kernel_type, param);
        corr_out = corr(kx, ky);

    return corr_out::Matrix{T} where T <: Real
end

function kxcorr(x::AbstractMatrix, y::AbstractMatrix, kernel_type::AbstractString)

    # --Test Polynomial
        if occursin("SIM", uppercase(kernel_type))
            error("Cannot do a similarity score for a kernel correlation, please select another kernel function")
        end 

    # --Same size
        (x, y) = same_length(x, y);

    # --Find Corr
        param = 10.0; 
        @inline kx = kernel(x, kernel_type, param);
        @inline ky = kernel(y, kernel_type, param);
        corr_out = corr(kx, ky);

    return corr_out::Matrix{T} where T <: Real
end

