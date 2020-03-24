module Analytical

using PyCall
using SpecialFunctions
using Roots

mutable struct parameters
	gam_neg::Int64
	gL::Int64
	gH::Int64
	alLow::Float64
	alTot::Float64
	theta_f::Float64
	theta_mid_neutral::Float64
	al::Float64
	be::Float64
	B::Float64
	pposL::Float64
	pposH::Float64
	N::Int64
	n::Int64
	Lf::Int64
	rho::Float64
	L_mid::Int64
	al2::Float64
	be2::Float64
	gF::Int64
	ABC::Bool
end
mutable struct popSizeParameters
	NN::Int64
	nn::Int64
end
adap = parameters(-83,10,500,0.2,0.2,0.001,0.001,0.184,0.000402,0.999,0.001,0.0,500,25,10^6,0.001,501,0.00415,0.00515625,500,false)
popSize = popSizeParameters(adap.N*2,adap.n*2)

export adap

################################
###### Solving parameters ######
################################

function Br(Lmax,theta)

	NN = popSize.NN
	gam_neg = adap.gam_neg
	rho = adap.rho
	t = -1.0*adap.gam_neg/(popSize.NN+0.0)
	u = theta/(2.0*popSize.NN)
	r = rho/(2.0*popSize.NN)
	#Lmax = 1000000

	return exp(-4*u*Lmax/(2*Lmax*r+t))
end

function set_Lf()
	B=0.999
	theta_f=5.417709306354578e-7
	i(L) = Br(L,theta_f)-B
	tmpLf  = find_zero(i,100)
	#tmpLf  = sc.fsolve(lambda L: Br(L,theta_f)-B,100)
	#Lf = convert(Int64,round(tmpLf[0]))
	Lf = convert(Int64,round(tmpLf))
	return Lf
end

function set_theta_f()
	i(theta) = Br(adap.Lf,theta)-adap.B
	#theta_f  = fsolve(lambda theta: Br(Lf,theta) - B,0.00001)
	theta_f  = find_zero(i,0.00001)
	adap.theta_f = theta_f
end

function alphaExpSimLow(pposL,pposH)
	return fixPosSim(adap.gL,0.5*pposL)/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH) + fixNegB(0.5*pposL+0.5*pposH))
end

function alphaExpSimTot(pposL,pposH)

	return (fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH))/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH)+fixNegB(0.5*pposL+0.5*pposH))
end

function solvEqns(params)
	pposL,pposH = params
	return (alphaExpSimTot(pposL,pposH)-adap.alTot,alphaExpSimLow(pposL,pposH)-adap.alLow)
end

function setPpos()

	sc = pyimport("scipy.optimize")
	pposL,pposH =  sc.fsolve(solvEqns,(0.001,0.001))
	#pposL = ppos[1]; pposH = ppos[2]
	#pposL = find_zero(alphaExpSimTot(adap.pposL,adap.pposH)-adap.alTot,0.001)
	#pposH =  find_zero(alphaExpSimLow(adap.pposL,adap.pposH)-adap.alLow,0.001)
	if pposL < 0.0
	 	pposL = 0.0
	end
	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
	if (pposH < 0.0 || pposH < 1e-15)
		pposH = 0.0
	end
	adap.pposL, adap.pposH = pposL, pposH
end

################################
######     Fixations      ######
################################

function pFix(gamma)
	s = gamma/(popSize.NN+0.0)
	pfix = (1.0-exp(-2.0*s))/(1.0-exp(-2.0*gamma))
	if s >= 0.1
		pfix = exp(-(1.0+s))
		lim = 0
		while(lim < 200)
			pfix = exp((1.0+s)*(pfix-1.0))
			lim +=1
		pfix = 1-pfix
		end
	end
	return pfix
end

function fixNeut()
	return 0.255*(1.0/(adap.B*popSize.NN))
end

function fixNegB(ppos::Float64)
	return 0.745*(1-ppos)*(2^(-adap.al))*(adap.B^(-adap.al))*(adap.be^adap.al)*(-SpecialFunctions.zeta(adap.al,1.0+adap.be/(2.0*adap.B))+SpecialFunctions.zeta(adap.al,0.5*(2-1.0/(adap.N*adap.B)+adap.be/adap.B)))
end

function fixPosSim(gamma,ppos)

	S = abs(adap.gam_neg/(1.0*popSize.NN))
	r = adap.rho/(2.0*popSize.NN)
	u = adap.theta_f/(2.0*popSize.NN)
	s = gamma/(popSize.NN*1.0)

	p0 = SpecialFunctions.polygamma(1,(s+S)/r)
	p1 = SpecialFunctions.polygamma(1,(r+adap.Lf*r+s+S)/r)
	CC = 1.0

	return 0.745*ppos*exp(-2.0*S*u*(p0-p1)*CC^2/r^2)*pFix(gamma)
end


################################
######    Polymorphism    ######
################################

function DiscSFSSelNegDown(ppos)

	result = adap.B*(adap.theta_mid_neutral)*0.745*(binomOp()'DiscSFSSelNeg(ppos))
	return result[1:lastindex(result)-1]
end

function DiscSFSSelNeg(ppos)

	beta = adap.be/(1.0*adap.B)
	NN2 = convert(Int64, round(popSize.NN*adap.B, digits=0))

	xa = [round(i/(NN2+0.0),digits=6) for i in 0:NN2]
	z(x,ppos=ppos) = (1.0-ppos)*(2.0^-adap.al)*(beta^adap.al)*(-SpecialFunctions.zeta(adap.al,x+beta/2.0) + SpecialFunctions.zeta(adap.al,(2+beta)/2.0))/((-1.0+x)*x)
	solveZ = xa .|> z
	return 1.0/(popSize.NN+0.0).*solveZ
end

function binomOp()

	sc = pyimport("scipy.stats")
	NN2 = convert(Int64, round(popSize.NN*adap.B, digits=0))
	samples = permutedims([i for j in 0:NN2, i in 1:popSize.NN+1])

	samplesFreqs = permutedims([round(j/(NN2+0.0),digits=6) for j in 0:NN2, i in 1:popSize.nn+1])

	#return samples,popSize.NN,samplesFreqs
	return sc.binom.pmf(samples,popSize.NN,samplesFreqs)
end

function cumulativeSfs(sfsTemp)

	out  = Array{Float64}(undef, length(sfsTemp))
	out[1] = sum(sfsTemp)

	for i in 2:length(sfsTemp)
		app = out[i-1]-sfsTemp[i-1]
		if app > 0.0
			out[i] = app
		else
			out[i] = (0.0)
		end
	end
	return out
end
end # module
