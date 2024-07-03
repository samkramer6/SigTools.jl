#=  kernelcorrtests.jl

    Test cases for the kernel_corr.jl functions

    Sam Kramer
    April 7th, 2024
=#

include("../src/SigTools.jl")
using .SigTools
using Test
using JET

#########################################################################################
#                                       Test                                            #
#########################################################################################

@testset "kernelcorrtests.jl" begin

    # --Setup
        x = vec(ones(Float64, 1, 1000));
        y = vec(ones(Float64, 1, 2000));
        a = ones(Float64, 100, 100);
        b = ones(Float64, 100, 200);
        c = ones(Float64, 50, 150);
        d = ones(Float64, 200, 200);
    
    # --KCC testing
        @test_nowarn corr_out = kxcorr(x, "Poly", 5);
        @test length( kxcorr(x, "Poly", 5)) == length(x);
        @test typeof( kxcorr(x, "Poly", 5)) == Vector{Float64}
        
        @test_nowarn kxcorr(x, y, "rbf", 10);
        @test typeof( kxcorr(x, y, "rbf", 1.1)) == Vector{Float64}
        @test_throws ErrorException("Cannot do a similarity score for a kernel correlation, please select another kernel function") kxcorr(x, y, "sim", 10)
        @test length( kxcorr(x, y, "rbf", 1) ) == length(larger(y, x)); 
        
        @test_nowarn kxcorr(a, b, "Poly", 1.0);
        @test_nowarn kxcorr(c, d, "poly", 10);
        @test size( kxcorr(c, d, "Poly", 1.0) ) == size(d);
        @test typeof( kxcorr(c, d, "RBF", 10) ) == Matrix{Float64}
        

end

#########################################################################################
#                                                                                       #
#########################################################################################

# --Type safety
    a = ones(Float64, 100, 100);
    b = ones(Float64, 150, 200);
    x = vec(ones(Float64, 1, 50000));

    display(methods(kxcorr))
    println(" ")

# --Timing
    println("===== Timing =====")               # It is much faster to specify a parameter
    @time kxcorr(a, b, "Poly", 3);
    @time kxcorr(a, b, "Poly");                 # Passing a parameter is much faster
    @time kxcorr(x, "RBF", 10);                 # Float64 much faster than Int64
    @time kxcorr(x, x, "RBF", 10);
