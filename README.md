# SigTools

This is Julia SigTools, a digital signal processing toolbox that was made to fill the gaps of other toolboxes. and streamline signal processing coding.

## Signal Generation

SigTools.jl comes with multiple different signal generation tools which can be used to create standard signals.

### Chirp functions (Swept Frequency Cosine Function)
- Able to create a swept frequency cosine wave function where the user is able to customize the frequency sweep rate as a function of time
- Comes in a form of logarithmic or linear sweeps

```julia
signal::Chirp = chirp(f_start::T where T<: Real,        # Start Frequency 
                      f_stop::T where T<: Real,         # End Frequency 
                      T::T where T<: Real,              # Period (s) 
                      cycles::Int64,                    # Number of cycles 
                      type::String,                     # F sweep type
                      fs::Int64,                        # Sample Rate 
                     )
```

where type `Chirp` is a custom type with fields of `signal` and `time` both of which are of type `Vector{Float64}`.

Passing `type::String = "Linear"` will create a linear frequency sweep pattern, where passing `type::String = "Log"` will generate an exponentially increasing frequency pattern. cycles defaults to 1. type `methods(chirp)` to see a list of comprehensive methods. 

### Heaviside step function

### PWM Signal

### Sigmoid

### Gaussian Pulse

## Noise Generation

SigTools.jl comes with different types of noise generation tools which come prepackaged or can be edited and customized. There are differnt types of colored and white noise genreators available. All of which are mean-zero and Gaussian.

### White Gaussian Noise
```julia
noise::Vector{Float64} = awgn(N::Int64,                           # Number of points in vector
                              A::T where {T <: Real},             # Amplitude of noise
                             ) 

noise::Vector{Float64} = white_noise(σ::T where {T<:Real},        # Amplitude Offset 
                                     N::Int64,                    # Number of points in vector 
                                     A::T where {T<:Real},        # Amplitude of noise
                                    )
```

These are the white noise generation functions provided by SigTools. `awgn()` is for mean-zero gaussian noise of a specified amplitude and length. `white_noise()` allows you to specify a mean as the amplitude offset as σ. Both provide a flat power spectral density characteristic of white noise.

### Colored Noise

Colored noise follows a white noise generation and filter process. This process is outlined Xu, Chang,\*An Easy Algorithm to Generate Colored Noise Sequences*, The Astronomical Journal, 2019. This function uses an Auto Regressive process by involving trucation of impulse response coefficients convolved with white noise. The filter order is 255.

## Utility Functions

## Signal Detections

## Window Functions
