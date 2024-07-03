#=  kerneltests.jl

    This is the test cases for kernels.jl

    Sam Kramer
=#

include("../src/SigTools.jl")
using .SigTools
using Test
using JET

#########################################################################################
#                                        Tests                                          #
#########################################################################################

@testset "kernels.jl" begin

    # --Setup
        x = collect(Float64, 1:.01:100);
        y = collect(Float64, 1:.01:100);
        a = ones(Float64, 100, 100);
        b = ones(Float64, 100, 200);
        c = ones(Float64, 200, 50);

    # --Linkernel tests
        @test_nowarn linkernel(x);
        @test typeof( linkernel(y) ) == Vector{Float64}
        @test length( linkernel(x) ) == length(x)

        @test_nowarn linkernel(x, y)
        @test typeof( linkernel(x, y) ) == Vector{Float64}
        @test length( linkernel(x, y) ) == length(y)

        @test_nowarn linkernel(x, y, "simscore")
        @test length( linkernel(x, y, "sim") ) == 1
        @test typeof( linkernel(x, y, "sim") ) == Float64
    
    # --Kernel tests
        @test_nowarn kernel(x, "RBF")
        @test_nowarn kernel(x, y, "Poly")
        @test_nowarn kernel(a, b, "RBF")
        @test_nowarn kernel(a, "RBF", 10)
        @test_nowarn kernel(a, c, "COSSIM", 1.0)
        @test_nowarn kernel(a, c, "Laplacian", 1)
        @test size( kernel(x, "Poly") ) == size(x)
        @test size( kernel(a, b, "Poly") ) == size(b)
        @test typeof( kernel(a, b, "RBF", 10) ) == Matrix{Float64}
end

#########################################################################################
#                                      Safety and Speed                                 #
#########################################################################################

# --Safety
    println("====== Saftey ======")
    display( methods(linkernel) )
    x = collect(Float64, 1:0.1:10000);
    display( @report_opt linkernel(x) )
    display( @report_opt linkernel(x, x) )
    display( @report_opt linkernel(x, x, "normsimscore") )
    
    display( methods(kernel) )
    display( @report_opt kernel(x, "RBF") )
    display( @report_opt kernel(x, x, "Poly") )
    display( @report_opt kernel(x, "RBF", 1.0) )
    display( @report_opt kernel(x, x, "Poly", 1) )
         
# --Timing
    println("====== Timing ======")
    @time linkernel(x)
    @time linkernel(x, x)
    @time linkernel(x, x, "abs")
    @time kernel(x, "RBF")
    @time kernel(x, x, "Poly")
    @time kernel(x, "Poly", 1.0)
    @time kernel(x, x, "RBF", 1.0)
