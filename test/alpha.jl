using Analytical
using Test


# Running Analytical approximation with default parameters
@testset "Analytical.jl" begin
    # @test changeParameters(-83,10,200,0.2,0.2,0.001,0.001,0.184,0.000402,0.4,0.001,0.0,500,25,10^6,501,0.001,0.0,0.0,5.0,false) == false
	Analytical.changeParameters(N=700,n=661,diploid=true,convoluteBinomial=true)
	Analytical.set_theta_f()
	Analytical.setPpos()
	@test adap.pposH == 0.0
	@test adap.theta_f == 5.417709306355064e-7

end
