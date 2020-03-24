using Analytical
using Test

@testset "Analytical.jl" begin
    @test Analytical.popSize.NN == 1000
    @test adap.gH == 200
end
