using Analytical
using Test

# Running Analytical approximation with default parameters
@testset "Analytical.jl" begin

	adap = Analytical.parameters(N=500,n=661,alLow=0.2)
	cnv = Analytical.binomialDict()
	Analytical.binomOp!(adap,cnv.bn)
	Analytical.setThetaF!(adap)
	Analytical.setPpos!(adap)

	@test adap.pposH == 0.0002641360146660052
	@test adap.thetaF == 1.6433217979109007e-6

end
