# Heaviside functions class testing files

include("../src/SigTools.jl")
using .SigTools
using Test
using JET

#########################################################################################
#                                       Tests                                           #
#########################################################################################

@testset "heaviside.jl" begin

  # --Setup
    (x, time) = stepfunction(1, 2);
    (a, ~) = pwmsignal(10, 50, 1.0);
    y = collect(-10:.01:10);

  # --Testing
    @test_nowarn stepfunction(1, 2);
    @test length(x) == length(time);
    @test maximum(x) == 1;
    @test minimum(x) == 0;
    @test x[1] == 0;
    @test typeof(time) == Vector{Float64}

  # --PWM signal testing
    @test_nowarn pwmsignal(10, 50, 1.0)
    @test maximum(a) == 1;
    @test minimum(a) == 0;
    @test a[1] == 0;
    @test a[end] == 0;

  # --Sigmoid testing
    @test_nowarn sigmoid(y, 1.0, 1.0, 1.0);
    @test_throws ErrorException("Cannot pass empty x vector") sigmoid([], 1.0, 1.0, 1.0);

  # --Gaussian testing
    @test maximum(gaussian(y, 10.0)) == 1; 
    @test_nowarn gaussian(y, 10.0);
    @test_throws ErrorException("Cannot pass empty x vector") gaussian([], 10.0);
end

##########################################################################################
#                                       Safety                                           #                             
##########################################################################################

x = collect(-10:0.01:10);
println("======== Safety Tests ========")
display(@report_opt pwmsignal(10, 50, 1.0))
display(@report_opt stepfunction(1, 2))
display(@report_opt sigmoid(x, 1.0, 1.0, 1.0))
display(@report_opt gaussian(x, 1.0))

println("======== Time testing ========")
@time pwmsignal(10, 50, 1.0)
@time stepfunction(1, 2)
@time sigmoid(x, 1.0, 1.0, 1.0)
@time gaussian(x, 1.0)

