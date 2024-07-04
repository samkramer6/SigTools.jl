# SigTools

This is Julia SigTools, a digital signal processing toolbox that was made to fill the gaps of other toolboxes. and streamline signal processing coding.

[![Build Status](https://github.com/samkramer6/SigTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/samkramer6/SigTools.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Signal Generation

SigTools.jl comes with multiple different signal generation tools which can be used to create standard signals.

### Chirp functions (Swept Frequency Cosine Function)
- Able to create a swept frequency cosine wave function where the user is able to customize the frequency sweep rate as a function of time
- Comes in a form of logarithmic or linear sweeps

```julia
signal::Chirp = chirp(f_start::T where T<: Real, 
                      f_stop::T where T<: Real, 
                      T::T where T<: Real, 
                      cycles::Int64, 
                      type::String,
                      fs::Int64,
                     )
```

where type `Chirp` is a custom type with fields of `signal` and `time` both of which are of type `Vector{Float64}`.

Passing `type::String = "Linear"` will create a linear frequency sweep pattern, where passing `type::String = "Log"` will generate an exponentially increasing frequency pattern. cycles defaults to 1. type `methods(chirp)` to see a list of comprehensive methods. 

### Heaviside step function

### PWM Signal

### Sigmoid

### Gaussian Pulse

## Noise Generation

## Utility Functions

## Signal Detections

## Window Functions
