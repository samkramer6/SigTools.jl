#=  heaviside.jl

    This file is going to describe the heaviside step functions for the SigTools.jl package

    Sam Kramer
    March 29th, 2024
=#

"""
SigTools.jl step()
This is the step function definition

Inputs:
        

Outputs:
        
    
Sam Kramer    
Methods:
"""
function stepfunction(start::Int64, stop::Int64)

  # --Start
  fs = 10000
  x = vec([0.0])
  if start != 0
    off_start = vec(zeros(1, round(fs * start)))
    append!(x, off_start)
  end

  # --On section
  on_start = vec(ones(1, fs .* (stop - start)))
  append!(x, on_start)

  # --Off section
  off_end = vec(zeros(1, round(fs * start)))
  append!(x, off_end)

  # --Time vector
  time = vec(collect(0:1/fs:(length(x)-1)./fs))

  return x::Vector{Float64}, time::Vector{Float64}
end

function pwmsignal(timelength::Int, dutycycle::Int, pulselen::Float64)

  # --Setup
  fs = 10000
  x = vec([0.0])
  time_on = timelength * (dutycycle / 100)
  num_cycles = floor(time_on ./ pulselen)
  time_off = timelength * (100 - dutycycle) / 100
  pulse_off = time_off / num_cycles

  # --Setup
  for i in collect(1:num_cycles)

    pulse = vec(ones(1, round(Int, fs .* pulselen)))
    append!(x, pulse)
    no_pulse = vec(zeros(1, round(Int, fs .* pulse_off)))
    append!(x, no_pulse)

  end

  # --Create time
  time = vec(collect(0:1/fs:(length(x)-1)./fs))

  return x::AbstractVector, time::AbstractVector
end

function sigmoid(x::AbstractVector, a::Float64, b::Float64, c::Float64)

  # --Find length to make sure its not empty x
  if length(x) === 0
    error("Cannot pass empty x vector")
  end

  # --calculate sigmoid function
  @fastmath y = a .* 1 ./ (1 .+ exp.(-b .* (x .- c)))

  return y::AbstractVector
end

function gaussian(x::AbstractVector, σ::Float64)

  # --Find length to make sure its not empty x
  if length(x) === 0
    error("Cannot pass empty x vector")
  end

  # --Calculate the gaussian function
  @fastmath A = 1 ./ (σ .* sqrt.(2 .* π))
  @fastmath y = exp.(-(x .^ 2) ./ (2 .* σ))

  return y::AbstractVector
end
