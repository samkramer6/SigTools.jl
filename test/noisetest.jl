include("../src/SigTools.jl")
using .SigTools
using Test

@testset "noisetest.jl" begin
    N = 2056
    A = 1

    @test_nowarn awgn(N, A)
    @test_nowarn white_noise(0, N, A)
    @test_nowarn red_noise(N, A)
    @test_nowarn deep_red_noise(N, A)
    @test_nowarn voss_noise(N, A)

    @test length(awgn(N, A)) == N
    @test length(white_noise(0, N, A)) == N
    @test length(red_noise(N, A)) == N
    @test length(deep_red_noise(N, A)) == N
    @test length(voss_noise(N, A)) == N

end
