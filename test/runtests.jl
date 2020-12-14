using Analytical
using Test

# Running Analytical approximation with default parameters
@testset "Analytical.jl" begin

    adap = Analytical.parameters(N=500,n=661,alLow=0.2)
    Analytical.binomOp!(adap)
    Analytical.setThetaF!(adap)
    Analytical.setPpos!(adap)

    @test adap.pposH == 0.00026413601466600506
    @test adap.thetaF == 1.6433217979109018e-6

end
