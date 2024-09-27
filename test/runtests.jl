include("../src/SigTools.jl")
using .SigTools
using Test
using StatsBase
using Pkg

Pkg.resolve()
Pkg.activate()
Pkg.instantiate()

running_tests_message = "
#####################################################################################
#                                    Running Tests                                  #
#####################################################################################

"
println(running_tests_message)

@testset "convolutions.jl" begin

  # --Setup
  a = vec(ones(Float64, 1, 1000))
  b = vec(ones(Float64, 1, 2000))
  c = ones(Int64, 100, 100)
  d = ones(Int64, 200, 200)
  e = ones(Float64, 200, 100)
  f = ones(Float64, 100, 200)

  # --larger() tests
  @test_nowarn a = larger(a, b)
  @test length(larger(a, b)) == length(b)

  # --same_length() tests
  @test_nowarn (a, b) = same_length(a, b)
  @test length(a) == length(b)
  @test_nowarn (c, d) = same_length(c, d)
  @test length(c) == length(d)
  @test_nowarn (d, e) = same_length(d, e)
  @test size(d, 1) == size(e, 1)
  @test size(e, 2) == size(d, 2)
  @test_nowarn (e, f) = same_length(e, f)
  @test size(e, 1) == size(f, 1)
  @test size(e, 2) == size(f, 2)

  # --Convolutions
  @test_nowarn convol(c, d)
  @test_nowarn convol(a, b)

end

@testset "correlations.jl" begin

  # --Setup tests
  x = vec(rand(Float64, 1, 1000))
  y = vec(rand(Int64, 1, 2000))
  e = [1 2; 3 4; 5 6]
  f = [1 2 3; 3 4 5; 6 7 8]

  # --1D correlation tests
  @test length(corr(x)) == length(x)
  @test length(corr(x, y)) == length(y)
  @test length(corr(y)) == length(y)
  @test length(corr(x, y, "Normalized")) == length(x)
  @test_throws ErrorException("Not a recognized type of correlation, options are: None, Normalized, Biased, and Unbiased.") corr(x, y, "ThrowError")

  # --2D correlation tests
  @test length(corr(e, f)) == length(f)
  @test size(corr(e, f), 1) == size(f, 1)
  @test size(corr(e, f), 2) == size(f, 2)

end

@testset "kernels.jl" begin

  # --Setup
  a = vec(ones(1, 1000))
  b = vec(rand(Float64, 1, 1000))
  c = vec(rand(Int64, 1, 2000))
  d = ones(Int64, 100, 100)
  e = ones(Float64, 100, 100)

  # --1D Linear Kernel functions
  @test_nowarn linkernel(a)
  @test length(linkernel(a)) == length(a)
  @test maximum(linkernel(a)) == 1
  @test linkernel(a, a, "sim") == sum(a)
  @test length(linkernel(a, b, "similarity")) == 1
  @test maximum(linkernel(a, a, "Random")) == 2
  @test length(linkernel(a, c, "NormSimScore")) == 1

  # --2D Linear Kernel Functions 
  @test size(linkernel(d)) == size(d)
  @test mean(linkernel(d)) == 1
  @test length(linkernel(d, d, "NormSimScore")) == 1
  @test size(linkernel(d, d)) == size(d)

  # --1D Nonlinear kernel functions
  @test_nowarn kernel(a, "rbf")
  @test length(kernel(b, c, "COSSIM")) == 1
  @test length(kernel(b, c, "poly")) == length(c)

  # --2D Nonlinear kernels
  @test length(kernel(d, "COSSIM")) == 1
  @test size(kernel(d, e, "Poly")) == size(e)
  @test_nowarn kernel(d, "RBF")

  # --Norms 
  @test length(euclidian_norm(a)) == 1
  @test euclidian_norm(a) == sqrt(sum(a))
  @test norm(a, 1) == sum(a)
  @test euclidian_norm(a) == norm(a, 2)
  @test infnorm(a) == 1
  @test infnorm(a) == 1.0

end

@testset "kernel_corr.jl" begin

  # --Setup
  a = vec(ones(Float64, 1, 1000))
  b = vec(zeros(Int64, 1, 10))
  c = ones(Float64, 10, 10)
  d = ones(Int64, 100, 100)

  # --kxcorr
  @test length(kxcorr(a, "poly", 10)) == length(a)
  @test_nowarn kxcorr(a, b, "poly", 10)
  @test_nowarn kxcorr(c, d, "Gaussian", 10)
  @test_throws ErrorException("Cannot do a similarity score for a kernel correlation, please select another kernel function") kxcorr(a, "Simscore", 10)
  @test_throws ErrorException("Cannot do a similarity score for a kernel correlation, please select another kernel function") kxcorr(a, b, "Simscore", 10)
end

@testset "logrange.jl" begin

  # --Test 1
  @test maximum(logrange(1.0, 1000.0, 30)) == 1000
  @test_nowarn logrange(1.0, 1000.0, 30)
  @test_nowarn logrange(1, 1000, 30)
  @test maximum(logrange(1, 1000, 1000)) == 1000
  @test length(logrange(1, 1000, 100)) == 100

end

@testset "chirps.jl" begin


end

@testset "heaviside.jl" begin

  # --Setup
  (x, time) = stepfunction(1, 2)
  (a, ~) = pwmsignal(10, 50, 1.0)

  # --Testing
  @test_nowarn stepfunction(1, 2)
  @test length(x) == length(time)
  @test maximum(x) == 1
  @test minimum(x) == 0
  @test x[1] == 0

  # --PWM signal testing
  @test_nowarn pwmsignal(10, 50, 1.0)
  @test maximum(a) == 1
  @test minimum(a) == 0
  @test a[1] == 0
  @test a[end] == 0
end
