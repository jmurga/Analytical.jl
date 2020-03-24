using Analytical
using Test


# Running Analytical approximation with default parameters
@testset "Analytical.jl" begin
    @test Analytical.popSize.NN == 1000
    @test adap.gH == 500
end
