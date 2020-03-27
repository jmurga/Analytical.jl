using Analytical
using Test


Running Analytical approximation with default parameters
@testset "Analytical.jl" begin
    # @test changeParameters(-83,10,200,0.2,0.2,0.001,0.001,0.184,0.000402,0.4,0.001,0.0,500,25,10^6,501,0.001,0.0,0.0,5.0,false) == false
	Analytical.changeParameters(gam_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=0.999,pposL=0.001,pposH=0,N=500,n=25,Lf=10^6,L_mid=501,rho=0.001,al2= 0.0415,be2=0.00515625,TE=5.0,ABC=false)
	Analytical.set_theta_f()
	Analytical.setPpos()
	@test adap.pposH == 0.0
	@test adap.theta_f == 5.417709306354578e-7

end
