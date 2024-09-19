include("../src/SigTools.jl")
using .SigTools
using Test

@testset "noisetest.jl" begin
  N = 1024
  A = 1

  @test_nowarn awgn(N, A)
  @test_nowarn white_noise(0, N, A)
  @test_nowarn red_noise(N, A)
  @test_nowarn pink_noise(N, A)

  #TODO: Complete testing of the types of noise and their PSDs are correct

end
