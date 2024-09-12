#=  walsh.jl

This file will be on the creation of Walsh Matrix for creating CDMA orthogonal codes

Sam Kramer
September 4th, 2024

=#

using LinearAlgebra
using Random

function hadamard(n::Int64)

    if n == 1
        H_n = [1 1; 1 -1]
    else
        H_n = [hadamard(n - 1) hadamard(n - 1); hadamard(n - 1) -hadamard(n - 1)]
    end

    return Matrix(H_n)::Matrix{Int64}
end

function sync_code(n_users)

    n_bits = ceil(Int64, log2(n_users))
    return hadamard(n_bits)[1:n_users, :]
end

function asynch_code(n_bits)
    rng = MersenneTwister(1234)
    return bitrand(rng, 2^n_bits)::BitVector
end
