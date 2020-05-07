################################
####### Define parameters ######
################################

import Parameters: @with_kw
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
	rho::Float64                 = 0.001
	TE::Float64                = 5.0

	NN::Int64 = 1000
	nn::Int64 = 500
	bn::Dict = Dict(bRange[i] => zeros(nn+1,NN) for i in 1:length(bRange))
end

adap = parameters()

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


################################
###### Solving parameters ######
################################

# π/π_0 ≈ ℯ^(-4μL/2rL+t)
function Br(Lmax::Int64,theta::Float64)

	γ_neg = adap.gam_neg
	ρ  	  = adap.rho
	t     = -1.0*γ_neg/(adap.NN+0.0)
	μ     = theta/(2.0*adap.NN)
	r     = ρ/(2.0*adap.NN)

	return ℯ^(-4*μ*Lmax/(2*Lmax*r+t))
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value ) in a locus.
function set_theta_f()

	i(θ)         = Br(adap.Lf,θ)-adap.B
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
