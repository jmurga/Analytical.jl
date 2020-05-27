using Analytical
using Test


# Running Analytical approximation with default parameters
@testset "Analytical.jl" begin
    Analytical.changeParameters()
    @test Analytical.adap.NN == 2000
    @test adap.gH == 500
end
