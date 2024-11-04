include("../src/SigTools.jl")
using Test
using .SigTools

@testset "windowtest.jl" begin
  N = 1024
  a = 2
  b = 0.5

  @test_nowarn rectangular_window(N)
  @test_nowarn triangle_window(N)
  @test_nowarn hann(N, a)
  @test_nowarn hamming(N)
  @test_nowarn blackman(N)
  @test_nowarn riesz(N)
  @test_nowarn riemann(N)
  @test_nowarn tukey(N, b)
  @test_nowarn poisson(N, b)

end
