module Analytical

using Parameters
using PyCall
using SpecialFunctions
using Roots

@with_kw mutable struct parameters
	gam_neg::Int64 = -83
	gL::Int64 = 10
	gH::Int64 = 500
	alLow::Float64 = 0.2
	alTot::Float64 = 0.2
	theta_f::Float64 = 1e-3
	theta_mid_neutral::Float64 = 1e-3
	al::Float64 = 0.184
	be::Float64 = 0.000402
	B::Float64 = 0.999
	pposL::Float64 = 0.001
	pposH::Float64 = 0
	N::Int64 = 500
	n::Int64 = 250
	Lf::Int64 = 10^6
	L_mid::Int64 = 501
	rho::Float64 = 0.001
	al2::Float64 =  0.0415
	be2::Float64 = 0.00515625
	ABC::Bool = false
end

@with_kw mutable struct popSizeParameters
	NN::Int64 = 1000
	nn::Int64 = 500
end

adap = parameters()
popSize = popSizeParameters(adap.N*2,adap.n*2)
export adap

################################
###### Solving parameters ######
################################

function Br(Lmax::Int64,theta::Float64)

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

function binomOp()

	sc = pyimport("scipy.stats")
	NN2 = convert(Int64, round(popSize.NN*adap.B, digits=0))

	samples =  [i for i in 0:popSize.nn, j in 0:NN2]
	samplesFreqs = [j for i in 0:popSize.nn, j in 0:NN2]
	samplesFreqs = samplesFreqs./NN2

	return sc.binom.pmf(samples,popSize.nn,samplesFreqs)
end

bn = binomOp()


################################
######     Fixations      ######
################################

function pFix(gamma::Int64)
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

function fixPosSim(gamma::Int64,ppos::Float64)

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

############Nuetral#############

function DiscSFSNeutDown()
	NN2 = convert(Int64,round(popSize.NN*adap.B))
	function neutralSfs(i)
		if i > 0 && i < NN2
			 return 1.0/(i)
		end
		return 0.0
	end
	x = [convert(Float64,i) for i in 0:NN2]
	solvedNeutralSfs = x .|> neutralSfs
	return adap.B*(adap.theta_mid_neutral)*0.255*((bn*solvedNeutralSfs))
end

############Positive############
# Variable gamma in function changed to gammaValue to avoid problem with exported SpecialFunctions.gamma

function DiscSFSSelPosDown(gammaValue, ppos)

	S = abs(adap.gam_neg/(1.0*popSize.NN))
	r = adap.rho/(2.0*popSize.NN)
	u = adap.theta_f/(2.0*popSize.NN)
	s = gammaValue/(popSize.NN*1.0)
	p0 = SpecialFunctions.polygamma(1,(s+S)/r)
	p1 = SpecialFunctions.polygamma(1,1.0+(r*adap.Lf+s+S)/r)
	red_plus = exp(-2.0*S*u*(p0-p1)/(r^2))

	# Solving sfs
	NN2 = convert(Int64,round(popSize.NN*adap.B,digits=0))
	x = [i for i in 0:NN2]
	x = x/(NN2+0.0)

	function positiveSfs(i,gammaCorrected=gammaValue*adap.B,ppos=ppos)
		if i > 0 && i < 1.0
			return ppos*0.5*(exp(2*gammaCorrected)*(1-exp(-2.0*gammaCorrected*(1.0-i)))/((exp(2*gammaCorrected)-1.0)*i*(1.0-i)))
		end
		return 0.0
	end

	solvedPositiveSfs = (1.0/(NN2+0.0)) * (x .|> positiveSfs)

	return (adap.theta_mid_neutral)*red_plus*0.745*(bn*solvedPositiveSfs
	)
end

######Slightly deleterious######
function DiscSFSSelNegDown(ppos)
	result = adap.B*(adap.theta_mid_neutral)*0.745*(bn*DiscSFSSelNeg(ppos))
	return result[1:lastindex(result)-1]
end

function DiscSFSSelNeg(ppos)

	beta = adap.be/(1.0*adap.B)
	NN2 = convert(Int64, round(popSize.NN*adap.B, digits=0))

	xa = [round(i/(NN2+0.0),digits=6) for i in 0:NN2]
	z(x,ppos=ppos) = (1.0-ppos)*(2.0^-adap.al)*(beta^adap.al)*(-SpecialFunctions.zeta(adap.al,x+beta/2.0) + SpecialFunctions.zeta(adap.al,(2+beta)/2.0))/((-1.0+x)*x)
	solveZ = xa .|> z
	if (solveZ[1] == Inf || isnan(solveZ[1]))
		solveZ[1] = 0.0
	end
	if (solveZ[lastindex(solveZ)] == Inf || isnan(solveZ[lastindex(solveZ)]))
		solveZ[lastindex(solveZ)] = 0.0
	end

	return 1.0/(popSize.NN+0.0).*solveZ
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

################################
###    Summary statistics    ###
################################

function alphaByFrequencies(gammaL,gammaH,pposL,pposH)

	ret = Array{Float64}(undef, popSize.nn - 1)
	sel = Array{Float64}(undef, popSize.nn - 1)
	# Fixation
	fN = adap.B*fixNeut()
	fNeg = adap.B*fixNegB(0.5*pposH+0.5*pposL)
	fPosL = fixPosSim(gammaL,0.5*pposL)
	fPosH = fixPosSim(gammaH,0.5*pposH)

	# Polymorphism
	neut = cumulativeSfs(DiscSFSNeutDown())
	selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
	selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
	selN = cumulativeSfs(DiscSFSSelNegDown(adap.pposH+adap.pposL))

	# for i in range(0,len(selH)):
	# 	sel.append((selH[i]+selL[i])+selN[i])
	# for i in range(0,self.nn-1):
	# 	ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))
	return (neut,selH,selL,selN,fN,fNeg,fPosL,fPosH)
end

# set_theta_f()
# setPpos()
#
# function test()
#
# 	for t in 1:10^6
# 		c,d,e,f,g,h,i,j = alphaByFrequencies(adap.gL,adap.gH,adap.pposH,adap.pposL)
#
# end
################################
#### Old functions equations ###
################################
#
# alphaExpSimLow(pposL,pposH) = fixPosSim(adap.gL,0.5*pposL)/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH) + fixNegB(0.5*pposL+0.5*pposH))
#
# alphaExpSimTot(pposL,pposH) = (fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH))/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH)+fixNegB(0.5*pposL+0.5*pposH))

end # module
