module SigTools

include("convolutions.jl")
export same_length
export smaller
export larger
export conv

include("package_greet.jl")
export package_greet
export package_list

include("correlations.jl")
export corr

include("kernels.jl")
export linkernel
export kernel
export euclidian_norm
export norm
export infnorm

include("kernel_corr.jl")
export kxcorr

include("logrange.jl")
export logrange

include("chirps.jl")
export chirp    # Chirp Function for generation of signal
export Chirp    # Data Type Struct of type Chirp

include("heaviside.jl")
export stepfunction
export pwmsignal
export sigmoid
export gaussian

include("windows.jl")
export rectangular_window
export triangle_window
export Hann
export Hamming
export Blackman
export Riesz
export Riemann
export Tukey
export Poisson

include("noise.jl")
export awgn

end
