#   psd.jl
using FFTW
using .SigTools

""" SigTools.jl psd.jl

This is a function that is used to compute the Power Spectral Density of a signal in the dB scale. Data is windowed using a Cos^2(x) window.

Inputs:
    data::Vector{T} where {T<:Real} == The data vector that is desired to find the PSD from

Outputs:
    psd_output::Vector{Float64} == The output of the PSD in dB scale

Example Usage:
            psd_output = psd(data)

Samuel Kramer

"""
function psd(data::Vector{T} where {T<:Real})

    data = data .* hann(length(data), 2)
    @fastmath psd_data = 2 .* abs.(fft(data) .* conj.(fft(data))) ./ length(data)

    return (20 .* log10.(psd_data))::Vector{Float64}
end

""" SigTools.jl psd.jl

This function gives a welch approximation of the PSD using the Welch 1967, "The Use of Fast Fourier Transform for the Estimation of Power Spectra" paper. It autodefaults to 8 separate sections with largest possible length.

Inputs:
    data::Vector{T} where {T<:Real} == The data vector you wish to find the PSD of.

Outputs:
    psd_out::Vector{Float64} == The welch approximation of the PSD of the data returned in dB.

Sam Kramer

"""
function welch(data::Vector{T} where {T<:Real})

    K = 8
    L = 2 * fld(length(data), K)
    overlap = fld(L, 2)
    solution = vec(zeros(Float64, 1, length(data)))
    @inline W = hamming(L)    # Window definition
    @fastmath U = (1 / L) .* sum(W .^ 2)

    for i in collect(range(0, K - 1))


        if (i * overlap + L) > length(data)
            @inbounds d_window = vec(data[(i * overlap + 1):length(data)])
            @inbounds d_window = d_window .* hamming(length(d_window))
        else
            @inbounds d_window = vec(data[(i * overlap + 1):(i * overlap + L)] .* W)
        end

        if i === 0
            @inbounds d_window =
                vcat(d_window, vec(zeros(1, length(data) - length(d_window))))
        else
            @inbounds d_window = vcat(vec(zeros(1, i * overlap)), d_window)
            @inbounds d_window =
                vcat(d_window, vec(zeros(1, length(data) - length(d_window))))
        end

        @fastmath solution = solution .+ (fft(d_window) .^ 2)

    end

    @fastmath solution = (L / (U * K)) .* abs.(solution)

    return (20 .* log10.(solution))::Vector{Float64}
end
