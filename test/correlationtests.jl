#=  correlationstests.jl

    This is the test cases for the correlations.jl file
    
    Sam Kramer
    April 7th, 2024
=#

include("../src/SigTools.jl")
using .SigTools
using Test
using JET

#########################################################################################
#                                      Tests                                            #
#########################################################################################

@testset "correlationstests.jl" begin

  # --Setup tests
      x = vec(rand(Float64, 1, 1000));
      y = vec(rand(Int64, 1, 2000));
      e = [1 2; 3 4; 5 6];
      f = [1 2 3; 3 4 5; 6 7 8];
      a = ones(Float64, 100, 100);
      b = ones(Float64, 150, 200);
      c = ones(Float64, 200, 100);

  # --1D correlation tests
      @test_nowarn corr(x);
      @test_nowarn corr(x, y);
      @test_nowarn corr(x, x);
      
      @test length( corr(x) ) == length(x);
      @test length( corr(x, y) ) == length(y);
      @test length( corr(y) ) == length(y);
      @test length( corr(x, y, "Normalized") ) == length(x);
      @test_throws ErrorException("Not a recognized type of correlation, options are: None, Normalized, Biased, and Unbiased.") corr(x, y, "ThrowError")
      @test typeof( corr(x, y) ) == Vector{Float64};

  # --2D correlation tests
      @test length(corr(e, f)) == length(f);
      @test size(corr(e, f), 1) == size(f, 1);
      @test size(corr(e, f), 2) == size(f, 2);
      @test_nowarn corr(a, b); 

end

#########################################################################################
#                                  Timing and Safety                                    #
#########################################################################################

# --Safety
    println("===== Methods =====")
    println(methods(corr))

# --Timing
    a = vec(ones(Float64, 1, 50000));
    b = vec(collect(1:10000));
    c = ones(Float64, 100, 100);
    d = ones(Float64, 200, 100);
    e = ones(Float64, 200, 200);

    @time corr(a, b);
    @time corr(a);
    @time corr(a, b, "Biased");
    @time corr(c, d);
    @time corr(d, e);
