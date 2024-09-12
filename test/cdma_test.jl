include("../src/SigTools.jl")
using .SigTools
using Test
using LinearAlgebra

@testset "cdma_test.jl" begin

  @test_nowarn hadamard(4)
  @test_nowarn sync_code(10)
  @test_nowarn asynch_code(8)

  @test dot(sync_code(4)[1, :], sync_code(4)[4, :]) === 0

  @test size(hadamard(4)) == (2^4, 2^4)
  @test typeof(hadamard(4)) == Matrix{Int64}
  @test maximum(hadamard(4)) == 1
  @test minimum(hadamard(4)) == -1

end
