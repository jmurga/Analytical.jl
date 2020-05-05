module Analytical

include("features.jl")
include("inferParams.jl")

using Parameters
using PyCall
using SpecialFunctions
using Distributions
using Roots
using HDF5

@with_kw mutable struct parameters
	gam_neg::Int64             = -83
	gL::Int64                  = 10
	gH::Int64                  = 500
	alLow::Float64             = 0.2
	alTot::Float64             = 0.2
	theta_f::Float64           = 1e-3
	theta_mid_neutral::Float64 = 1e-3
	al::Float64                = 0.184
	be::Float64                = 0.000402
	B::Float64                 = 0.999
	bRange::Array{Float64,1}   = [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.999]
	pposL::Float64             = 0.001
	pposH::Float64             = 0
	N::Int64                   = 500
	n::Int64                   = 250
	Lf::Int64                  = 10^6
	rho::Float64               = 0.001
	TE::Float64                = 5.0

	NN::Int64 = 1000
	nn::Int64 = 500
	bn::Dict = Dict(bRange[i] => zeros(nn+1,NN) for i in 1:length(bRange))
end

adap = parameters()

export adap
################################
###### Solving parameters ######
################################

function changeParameters(;gam_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,bRange=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.999],B=0.999,pposL=0.001,pposH=0,N=500,n=25,Lf=10^6,rho=0.001,TE=5.0,convoluteBinomial=true)

	adap.gam_neg           = gam_neg
	adap.gL                = gL
	adap.gH                = gH
	adap.alLow             = alLow
	adap.alTot             = alTot
	adap.theta_f           = theta_f
	adap.theta_mid_neutral = theta_mid_neutral
	adap.al                = al
	adap.be                = be
	adap.bRange            = bRange
	adap.B                 = B
	adap.pposL             = pposL
	adap.pposH             = pposH
	adap.N                 = N
	adap.n                 = n
	adap.Lf                = Lf
	adap.rho               = rho
	adap.TE                = TE

	adap.NN = N*2
	adap.nn = n*2

	if convoluteBinomial == true
		adap.bn = Dict(bRange[i] => binomOp(bRange[i]) for i in 1:length(bRange))
	end
end

function Br(Lmax::Int64,theta::Float64)

	gam_neg = adap.gam_neg
	rho     = adap.rho
	t       = -1.0*adap.gam_neg/(adap.NN+0.0)
	u       = theta/(2.0*adap.NN)
	r       = rho/(2.0*adap.NN)

	return exp(-4*u*Lmax/(2*Lmax*r+t))
end

function set_Lf()

	i(L)   = Br(L,theta_f,adap,adap)-B
	tmpLf  = Roots.find_zero(i,100)
	Lf     = convert(Int64,round(tmpLf))

	return Lf
end

function set_theta_f()

	i(theta)     = Br(adap.Lf,theta)-adap.B
	theta_f      = Roots.find_zero(i,0.00001)
	adap.theta_f = theta_f
end

function alphaExpSimLow(pposL::Float64,pposH::Float64)
	return fixPosSim(adap.gL,0.5*pposL)/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH) + fixNegB(0.5*pposL+0.5*pposH))
end

function alphaExpSimTot(pposL::Float64,pposH::Float64)
	return (fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH))/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH)+fixNegB(0.5*pposL+0.5*pposH))
end

function solvEqns(params)

	pposL,pposH = params
	return (alphaExpSimTot(pposL,pposH)-adap.alTot,alphaExpSimLow(pposL,pposH)-adap.alLow)
end

function setPpos()
 	sc          = pyimport("scipy.optimize")
	# pposL,pposH = sc.fsolve(solvEqns,(0.001,0.001))
	pposL,pposH = sc.fsolve(solvEqns,(0.0,0.0))

	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
	# function f!(F,x)
	# 	F[1] = alphaExpSimTot(x[1],x[2])-adap.alTot
	# 	F[2] = alphaExpSimLow(x[1],x[2])-adap.alLow
	# end
	#
	# pposL,pposH = nlsolve(f!,[0.0; 0.0]).zero

	if pposL < 0.0
	 	pposL = 0.0
	end
	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
	if (pposH < 0.0 || pposH < 9e-15)
		pposH = 0.0
	end
	adap.pposL,adap.pposH = pposL, pposH
end

function binomOp(B)

    NN2          = convert(Int64, round(adap.NN*B, digits=0))
    samples      =  [i for i in 0:adap.nn]
    samplesFreqs = [j for j in 0:NN2]
    samplesFreqs = permutedims(samplesFreqs/NN2)

    f(x) = Distributions.Binomial(adap.nn,x)
    z    = samplesFreqs .|> f

    # out  = Array{Float64}(undef,(nn+1,NN))
    out  = Array{Float64}(undef,(adap.nn+1,adap.NN))
    out  = Distributions.pdf.(z,samples)
    return out
end

################################
######     Fixations      ######
################################

function pFix(gamma::Int64)
	s    = gamma/(adap.NN+0.0)
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
	return 0.255*(1.0/(adap.B*adap.NN))
end

function fixNegB(ppos::Float64)
	return 0.745*(1-ppos)*(2^(-adap.al))*(adap.B^(-adap.al))*(adap.be^adap.al)*(-SpecialFunctions.zeta(adap.al,1.0+adap.be/(2.0*adap.B))+SpecialFunctions.zeta(adap.al,0.5*(2-1.0/(adap.N*adap.B)+adap.be/adap.B)))
end

function fixPosSim(gamma::Int64,ppos::Float64)

	S  = abs(adap.gam_neg/(1.0*adap.NN))
	r  = adap.rho/(2.0*adap.NN)
	u  = adap.theta_f/(2.0*adap.NN)
	s  = gamma/(adap.NN*1.0)

	p0 = SpecialFunctions.polygamma(1,(s+S)/r)
	p1 = SpecialFunctions.polygamma(1,(r+adap.Lf*r+s+S)/r)
	CC = 1.0

	return 0.745*ppos*exp(-2.0*S*u*(p0-p1)*CC^2/r^2)*pFix(gamma)
end

################################
######    Polymorphism    ######
################################

############Neutral#############

function DiscSFSNeutDown()

	NN2 = convert(Int64,round(adap.NN*adap.B))

	function neutralSfs(i)
		if i > 0 && i < NN2
			 return 1.0/(i)
		end
		return 0.0
	end

	x                = [convert(Float64,i) for i in 0:NN2]
	solvedNeutralSfs = x .|> neutralSfs
	out              = adap.B*(adap.theta_mid_neutral)*0.255*(adap.bn[adap.B]*solvedNeutralSfs)

	return 	out[2:lastindex(out)-1]
end

############Positive############
# Variable gamma in function changed to gammaValue to avoid problem with exported SpecialFunctions.gamma

function DiscSFSSelPosDown(gammaValue::Int64,ppos::Float64)

	if ppos == 0.0
		out = zeros(Float64,adap.nn + 1)
	else
		S        = abs(adap.gam_neg/(1.0*adap.NN))
		r        = adap.rho/(2.0*adap.NN)
		u        = adap.theta_f/(2.0*adap.NN)
		s        = gammaValue/(adap.NN*1.0)
		p0       = SpecialFunctions.polygamma(1,(s+S)/r)
		p1       = SpecialFunctions.polygamma(1,1.0+(r*adap.Lf+s+S)/r)
		red_plus = exp(-2.0*S*u*(p0-p1)/(r^2))

		# Solving sfs
		NN2 = convert(Int64,round(adap.NN*adap.B,digits=0))
		xa  = [i for i in 0:NN2]
		xa  = xa/(NN2+0.0)

		function positiveSfs(i,gammaCorrected=gammaValue*adap.B,ppos=ppos)
			if i > 0 && i < 1.0
				return ppos*0.5*(exp(2*gammaCorrected)*(1-exp(-2.0*gammaCorrected*(1.0-i)))/((exp(2*gammaCorrected)-1.0)*i*(1.0-i)))
			end
			return 0.0
		end

		solvedPositiveSfs = (1.0/(NN2+0.0)) * (xa .|> positiveSfs)
		out               = (adap.theta_mid_neutral)*red_plus*0.745*(adap.bn[adap.B]*solvedPositiveSfs)
	end

	return out[2:lastindex(out)-1]
end

######Slightly deleterious######
function DiscSFSSelNegDown(ppos::Float64)
	out = adap.B*(adap.theta_mid_neutral)*0.745*(adap.bn[adap.B]*DiscSFSSelNeg(ppos))
	return out[2:lastindex(out)-1]
end

function DiscSFSSelNeg(ppos::Float64)

	beta     = adap.be/(1.0*adap.B)
	NN2      = convert(Int64, round(adap.NN*adap.B, digits=0))
	xa       = [round(i/(NN2+0.0),digits=6) for i in 0:NN2]
	z(x,ppos=ppos) = (1.0-ppos)*(2.0^-adap.al)*(beta^adap.al)*(-SpecialFunctions.zeta(adap.al,x+beta/2.0) + SpecialFunctions.zeta(adap.al,(2+beta)/2.0))/((-1.0+x)*x)
	solveZ   = xa .|> z

	if (solveZ[1] == Inf || isnan(solveZ[1]))
		solveZ[1] = 0.0
	end
	if (solveZ[lastindex(solveZ)] == Inf || isnan(solveZ[lastindex(solveZ)]))
		solveZ[lastindex(solveZ)] = 0.0
	end

	return 1.0/(NN2+0.0).*solveZ
end

function cumulativeSfs(sfsTemp)

	out    = Array{Float64}(undef, length(sfsTemp)+1)
	out[1] = sum(sfsTemp)

	for i in 2:length(sfsTemp)+1
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
function poissonFixation(;observedValue, λds, λdn)

	poissonS  = (λds/(λds + λdn) .* observedValue) .|> Poisson
	poissonD  = (λdn/(λds + λdn) .* observedValue) .|> Poisson

	sampledDs = rand.(poissonS,1)
	sampledDn = rand.(poissonD,1)

	return(reduce(vcat,sampledDs),reduce(vcat,sampledDn))
end

function poissonPolymorphism(;observedValues, λps, λpn)

    psPois(x,y=λps,z=λpn) = rand.((y./(y .+ z) .* x) .|> Poisson,1)
    pnPois(x,y=λps,z=λpn) = rand.((z./(y .+ z) .* x) .|> Poisson,1)

    sampledPs = observedValues .|> psPois
    sampledPn = observedValues .|> pnPois

    return (reduce(hcat,sampledPs),reduce(hcat,sampledPn))
end

function alphaByFrequencies(gammaL::Int64,gammaH::Int64,pposL::Float64,pposH::Float64,observedData::Array,nopos::String)

	P = observedData[:,1][lastindex(observedData[:,1])]
	D = observedData[:,2][lastindex(observedData[:,2])]

	if nopos == "pos"
		# Fixation
		fN    = adap.B*fixNeut()
		fNeg  = adap.B*fixNegB(0.5*pposH+0.5*pposL)
		fPosL = fixPosSim(gammaL,0.5*pposL)
		fPosH = fixPosSim(gammaH,0.5*pposH)

		# Polymorphism
		neut = cumulativeSfs(DiscSFSNeutDown())
		selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
		selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))

		# Outputs
		sel = (selH+selL)+selN

		ret = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)
		ret = ret[1:lastindex(ret)-1]

		# ret = Array{Float64}(undef, adap.nn-1)
		# sel = Array{Float64}(undef, adap.nn - 1)
		# for i in 1:length(ret)
		# 	sel[i] = (selH[i]+selL[i])+selN[i]
		# 	ret[i] = float(1.0 - (fN/(fPosL + fPosH +  fNeg + 0.0)) * sel[i]/neut[i])
		# end

		return (ret)
	elseif nopos == "nopos"

		# Fixation
		fN    = adap.B*fixNeut()*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN
		fNeg  = adap.B*fixNegB(0.5*pposH+0.5*pposL)*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN
		fPosL = fixPosSim(gammaL,0.5*pposL)*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN
		fPosH = fixPosSim(gammaH,0.5*pposH)*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN

		ds = fN
		dn = fNeg + fPosL + fPosH

		# Polymorphism
		neut = cumulativeSfs(DiscSFSNeutDown())
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))
		sel = selN

		ps = sum(neut) ./ sum(sel+neut)
		pn = sum(sel) ./ sum(sel+neut)

		# Outputs
		expectedDs, expectedDn = poissonSampling(observedValue=D,lambdaS=ds,lambdaN=dn)
		expectedPs, expectedPn = poissonSampling(observedValue=P,lambdaS=ps,lambdaN=pn)

		ret = 1 .- (fN/(fPosL + fPosH+  fNeg+0.0)) .* (sel./neut)
		ret = ret[1:lastindex(ret)-1]

		expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,ret[lastindex(ret)])
		# expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,fill(adap.B,size(expectedPn)[1]),fill(ret[lastindex(ret)],size(expectedPn)[1]))

		return (ret,expectedValues)
	else

		## Accounting for positive alleles segregating due to linkage
		# Fixation
		fN     = adap.B*fixNeut()
		fNeg   = adap.B*fixNegB(0.5*pposH+0.5*pposL)
		fPosL  = fixPosSim(gammaL,0.5*pposL)
		fPosH  = fixPosSim(gammaH,0.5*pposH)

		ds = fN
		dn = fNeg + fPosL + fPosH

		## Polymorphism
		neut = cumulativeSfs(DiscSFSNeutDown())
		selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
		selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))

		sel = (selH+selL)+selN
		replace!(sel,NaN=>0)
		ps = sum(neut) ./ sum(sel+neut)
		pn = sum(sel) ./ sum(sel+neut)

		## Outputs
		expectedDs, expectedDn = poissonSampling(observedValue=D,lambdaS=ds,lambdaN=dn)
		expectedPs, expectedPn = poissonSampling(observedValue=P,lambdaS=ps,lambdaN=pn)

		ret = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)
		ret = ret[1:lastindex(ret)-1]

		expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,ret[lastindex(ret)])

		# Accounting only for neutral and deleterious alleles segregating
		## Fixation
		fN_nopos     = fN*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fNeg_nopos   = fNeg*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fPosL_nopos  = fPosL*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fPosH_nopos  = fPosH*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN

		ds_nopos = fN_nopos
		dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

		## Polymorphism
		sel_nopos = selN
		ps_nopos = sum(neut) ./ sum(sel_nopos+neut)
		pn_nopos = sum(sel_nopos) ./ sum(sel_nopos+neut)

		## Outputs
		expectedDs_nopos, expectedDn_nopos = poissonSampling(observedValue=D,lambdaS=ds_nopos,lambdaN=dn_nopos)
		expectedPs_nopos, expectedPn_nopos = poissonSampling(observedValue=P,lambdaS=ps_nopos,lambdaN=pn_nopos)

		ret_nopos = 1 .- (fN_nopos/(fNeg_nopos + fPosH_nopos+  fNeg_nopos+0.0)) .* (sel_nopos./neut)
		ret_nopos = ret_nopos[1:lastindex(ret_nopos)-1]

		expectedValues = hcat(expectedValues,ret_nopos[lastindex(ret_nopos)])

		return (ret,ret_nopos,expectedValues)
	end
end

function summaryStatistics(fileName,alpha,expectedValues)


	h5open(fileName, "cw") do file
		tmp = string(rand(Int64))
	    write(file, "tmp"*tmp*"/alpha", alpha)
	    write(file, "tmp"*tmp*"/expectedValues", expectedValues)
	end
end

end # module


################################
######## Old functions  ########
################################
# function changeParameters(;gam_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,bRange=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.999],B=0.95,pposL=0.001,pposH=0,N=500,n=25,Lf=10^6,L_mid=501,rho=0.001,al2= 0.0415,be2=0.00515625,TE=5.0,ABC=false,convoluteBinomial=true)
#
# 	adap.gam_neg           = gam_neg
# 	adap.gL                = gL
# 	adap.gH                = gH
# 	adap.alLow             = alLow
# 	adap.alTot             = alTot
# 	adap.theta_f           = theta_f
# 	adap.theta_mid_neutral = theta_mid_neutral
# 	adap.al                = al
# 	adap.be                = be
# 	adap.bRange            = bRange
# 	adap.B                 = B
# 	adap.pposL             = pposL
# 	adap.pposH             = pposH
# 	adap.N                 = N
# 	adap.n                 = n
# 	adap.Lf                = Lf
# 	adap.L_mid             = L_mid
# 	adap.rho               = rho
# 	adap.al2               = al2
# 	adap.be2               = be2
# 	adap.TE                = TE
# 	adap.ABC               = ABC
#
# 	adap.NN = N*2
# 	adap.nn = n*2
#
# 	if convoluteBinomial == true
# 		adap.bn = Dict(bRange[i] => binomOp(bRange[i]) for i in 1:length(bRange))
# 	end
# end
# function setPpos_nlsolve()
#
# 	function f!(F,x)
# 		F[1] = alphaExpSimTot(x[1],x[2])-adap.alTot
# 		F[2] = alphaExpSimLow(x[1],x[2])-adap.alLow
# 	end
#
# 	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
# 	pposL,pposH = NLsolve.nlsolve(f!,[ 0.0; 0.0]).zero
#
# 	if pposL < 0.0
# 	 	pposL = 0.0
# 	end
# 	if (pposH < 0.0 || pposH < 9e-15)
# 		pposH = 0.0
# 	end
# 	adap.pposL,adap.pposH = pposL, pposH
# end

# function summaryStatistics(fileName,simulationName,alphaPos,alphaNopos)
#
# 	if isfile(fileName)
# 		jldopen(fileName, "a+") do file
# 			group = JLD2.Group(file, simulationName*string(rand(Int64)))
# 			group["alphaPos"] = alphaPos
# 			group["alphaNopos"] = alphaNopos
# 		end
# 	else
# 		jldopen(fileName, "w") do file
# 			group = JLD2.Group(file, simulationName*string(rand(Int64)))
# 			group["alphaPos"] = alphaPos
# 			group["alphaNopos"] = alphaNopos
# 		end
# 	end
# end

# function alphaByFrequenciesNoPositive(gammaL,gammaH,pposL,pposH)
#
# 	fN = adap.B*fixNeut()*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
# 	fNeg = adap.B*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN*fixNegB(0.5*pposH+0.5*pposL)
# 	fPosL = fixPosSim(gammaL,0.5*pposL)*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
# 	fPosH = fixPosSim(gammaH,0.5*pposH)*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
#
# 	neut = cumulativeSfs(DiscSFSNeutDown())
# 	selN = cumulativeSfs(DiscSFSSelNegDown(pposL+pposH))
#
# 	# Outputs
# 	ret = Array{Float64}(undef, adap.nn - 1)
# 	sel = Array{Float64}(undef, adap.nn - 1)
# 	for i in 1:length(ret)
# 		sel[i] = (selH[i]+selL[i])+selN[i]
# 		ret[i] = float(1.0 - (fNNopos/(fPosLNopos + fPosHNopos+  fNegNopos+0.0))* selNopos[i]/neutNopos[i])
# 	end
#
# 	return ret
# end

# set_theta_f()
# setPpos()

# adap.bn = Analytical.binomOp()
# function test()
# 		Analytical.set_theta_f()
# 		theta_f = adap.theta_f
# 		adap.B = 0.999
# 		Analytical.set_theta_f()
# 		Analytical.setPpos()
# 		adap.theta_f = theta_f
# 		adap.B = 0.4
# 		c,d,e,f,g,h,i,j = Analytical.alphaByFrequencies(adap.gL,adap.gH,adap.pposH,adap.pposL)
# end


#
# alphaExpSimLow(pposL,pposH) = fixPosSim(adap.gL,0.5*pposL)/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH) + fixNegB(0.5*pposL+0.5*pposH))
#
# alphaExpSimTot(pposL,pposH) = (fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH))/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH)+fixNegB(0.5*pposL+0.5*pposH))

# function binomOp()
#
# 	sc = pyimport("scipy.stats")
# 	NN2 = convert(Int64, round(adap.NN*adap.B, digits=0))
#
# 	samples =  [i for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = [j for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = samplesFreqs./NN2
#
# 	return sc.binom.pmf(samples,adap.nn,samplesFreqs)
# 	NN2          = convert(Int64, round(adap.NN*adap.B, digits=0))
# 	samples      =  [i for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = [j for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = samplesFreqs/NN2
#
# 	f(x) = Distributions.Binomial(adap.nn,x)
# 	z    = samplesFreqs .|> f
#
# 	out  = Array{Float64}(undef,(adap.nn+1,adap.NN))
# 	out  = Distributions.pdf.(z,samples)
#
# 	return out
# end
