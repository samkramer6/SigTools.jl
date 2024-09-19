# using .SigTools
using FFTW
using Test

########################################################################
#                         1D cross-correlation                         #
########################################################################

function corr(x::Vector{Float64})

    f_x = fft(x)

    @fastmath corr = real.(ifft(f_x .* conj.(f_x)))

    return corr::Vector{Float64}
end

function corr(x::Vector{T} where {T<:Real})

    f_x = fft(x)

    @fastmath corr = real.(ifft(f_x .* conj.(f_x)))

    return corr::AbstractVector
end

function corr(x::Vector{T} where {T<:Real}, type::AbstractString)

    f_x = fft(x)

    if uppercase(type) == "NONE"
        prefix = 1

    elseif uppercase(type) == "BIASED"

        prefix = 1 / (length(x))

    elseif uppercase(type) == "UNBIASED"

        m = collect(1:length(x))
        @inbounds prefix = 1 ./ (length(x) .- m)

    elseif uppercase(type) == "NORMALIZED"

        @fastmath rxx = real(ifft(f_x .* conj(f_x)))

        prefix = 1 ./ sqrt(rxx[1] * rxx[1])

    else

        error("Not a recognized type of correlation, options are: None, Normalized, Biased, and Unbiased.")

    end

    corr = prefix .* real.(ifft(f_x .* conj.(f_x)))

    return corr::Vector{Float64}
end

function corr(x::Vector{Float64}, y::Vector{Float64})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @fastmath corr = real.(ifft(f_x .* conj.(f_y)))

    return corr::Vector{Float64}
end

function corr(x::Vector{T} where {T<:Real}, y::Vector{T} where {T<:Real})

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    @fastmath corr = real.(ifft(f_x .* conj.(f_y)))

    return corr::AbstractVector
end

function corr(x::Vector{T} where {T<:Real}, y::Vector{T} where {T<:Real}, type::AbstractString)

    (x, y) = same_length(x, y)

    f_x = fft(x)
    f_y = fft(y)

    if uppercase(type) == "NONE"

        prefix = 1

    elseif uppercase(type) == "BIASED"

        prefix = 1 / (length(x))

    elseif uppercase(type) == "UNBIASED"

        m = collect(1:length(x))
        @inbounds prefix = 1 ./ (length(x) .- m)

    elseif uppercase(type) == "NORMALIZED"

        @fastmath rxx = real(ifft(f_x .* conj(f_x)))
        @fastmath ryy = real(ifft(f_y .* conj(f_y)))

        prefix = 1 ./ sqrt(rxx[1] * ryy[1])

    else

        error("Not a recognized type of correlation, options are: None, Normalized, Biased, and Unbiased.")

    end

    @fastmath corr = prefix .* real.(ifft(f_x .* conj.(f_y)))

    return corr::Vector{T} where {T<:Real}
end

########################################################################
#                         2D cross-correlation                         #
########################################################################

function corr(x::Matrix{Float64}, y::Matrix{Float64})

    # --Get same Length
    (x, y) = same_length(x, y)

    # --Find Gaussian Prefix
    m = ones(size(x, 1), size(x, 2))
    @fastmath prefix = 1 ./ exp.(-1 .* (m .^ 2) ./ (2 * length(x)))

    # --FFTs of data
    f_x = fft(x)
    f_y = fft(y)

    # --Find Correlation
    @fastmath corr = prefix .* real.(ifft(f_x .* conj.(f_y)))

    return corr::Matrix{Float64}
end

function corr(x::Matrix{T} where {T<:Real}, y::Matrix{T} where {T<:Real})

    # --Get same Length
    (x, y) = same_length(x, y)

    # --Find Gaussian Prefix
    m = ones(size(x, 1), size(x, 2))
    @fastmath prefix = 1 ./ exp.(-1 .* (m .^ 2) ./ (2 * length(x)))

    # --Find FFTs
    f_x = fft(x)
    f_y = fft(y)

    # --Find Correlation
    @fastmath corr = prefix .* real.(ifft(f_x .* conj.(f_y)))

    return corr::Matrix{T} where {T<:Real}
end

