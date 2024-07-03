#= convolutionstests.jl

    Test cases for functions within convolutions.jl.

    Test cases all pass, and working on optimization of code. 
    Notes:
        - Int64 is unstable for same_length and for conv(), it uses large amounts of memory
        - fft() package is unstable and might be solved by more strict type declarations
        - fft() already preplans, so it is a waste to preplan fft
=#

include("../src/SigTools.jl")
using .SigTools
using Test
using StatsBase
using JET

#########################################################################################
#                                     Test Cases                                        #
#########################################################################################

@testset "convolutions.jl" begin

  # --Setup
    a = vec(ones(Float64, 1, 1000));
    b = vec(ones(Float64, 1, 2000));
    c = ones(Int64, 100, 100);
    d = ones(Int64, 200, 200);
    e = ones(Float64, 200, 100);
    f = ones(Float64, 100, 200);

  # --larger() tests
    @test_nowarn a = larger(a, b);
    @test length(larger(a, b)) == length(b);

  # --same_length() tests
    @test_nowarn (a, b) = same_length(a, b);
    @test length(a) == length(b);
    @test_nowarn (c, d) = same_length(c, d);
    @test typeof(c) == Matrix{Float64};
    @test length(c) == length(d);
    @test_nowarn (d, e) = same_length(d, e);
    @test size(d, 1) == size(e, 1);
    @test size(e, 2) == size(d, 2);
    @test_nowarn (e, f) = same_length(e, f);
    @test size(e, 1) == size(f, 1);
    @test size(e, 2) == size(f, 2);

  # --Convolutions
    @test_nowarn conv(c, d);
    @test_nowarn conv(a, b);

end

#########################################################################################
#                                         Safety Tests                                  #
#########################################################################################

# --Setup
    a = vec(ones(Float64, 1, 1000));
    b = vec(ones(Float64, 1, 2000));
    c = ones(Int64, 100, 100);
    d = ones(Int64, 200, 200);
    e = ones(Float64, 100, 100);
    f = ones(Float64, 200, 200);

# --Safety Tests
    println("========= Safety Tests =========")
    display(@report_opt same_length(a, b))
    display(@report_opt same_length(e, f)) #<-- Int64 is unstable for these calculations
    # display(@report_opt conv(a, b)) <-- Inherently fft() is type unstable so will throw
    display(@report_opt larger(a, b))
    display(@report_opt smaller(c, d, 1))

# --Timing Tests
    println("========= Timing Tests ==========")
    @time same_length(a, b)
    @time same_length(e, f)
    @time conv(a, b)
    @time conv(c, d)
    @time conv(e, f)
