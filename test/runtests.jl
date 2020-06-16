using Analytical
using Test

# Running Analytical approximation with default parameters
@testset "Analytical.jl" begin

	adap = Analytical.parameters(N=500,n=661,alLow=0.2)
    Analytical.binomOp(adap)
	Analytical.set_theta_f(adap)
    Analytical.setPpos(adap)
    
    @test adap.pposH == 0.0
	@test adap.theta_f == 5.417709306355064e-7

end
