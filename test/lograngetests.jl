include("../src/SigTools.jl")
using .SigTools
using Test
using JET

@testset "lograngetests.jl" begin

    # --Tests
    @test_nowarn logrange(1, 100, 1000)
    @test typeof(logrange(1, 100, 1000)) == Vector{Float64}
    @test_nowarn logrange(1.5, 100.55, 1000)
    @test typeof(logrange(1.0, 100.00, 1000)) == Vector{Float64}
end

# --Testing Stability
display(methods(logrange))
display(@report_opt logrange(1, 100, 1000))
display(@report_opt logrange(1.0, 100.00, 1000))

# --Timing tests
@time logrange(1, 1000, 100000);
@time logrange(1.0, 1000.0, 100000);

