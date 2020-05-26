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
	rho::Float64               = 0.001
	TE::Float64                = 5.0
	diploid					   = false

	NN::Int64 = 1000
	nn::Int64 = 500

	bn::Dict = Dict(bRange[i] => zeros(nn+1,NN) for i in 1:length(bRange))
end

"""
	adap( 
		gam_neg::Int64,
		gL::Int64,	
		gH::Int64,
		alLow::Float64,
		alTot::Float64,
		theta_f::Float64,
		theta_mid_neutral::Float64,
		al::Float64,
		be::Float64,
		B::Float64,
		bRange::Array{Float64,1}
		pposL::Float64,
		pposH::Float64,
		N::Int64,
		n::Int64,
		Lf::Int64,
		rho::Float64,
		TE::Float64,
	)

Mutable structure containing the variables required to solve the analytical approach. All the functions are solve using the internal values of the structure. For this reason, *adap* is the only exported variable. Adap should be change before the perform the analytical approach, in other case, ```\$\\alpha_{(x)}\$``` will be solve with the default values.

# Parameters
 - `gam_neg::Int64`: 
 - `gL::Int64`: 
 - `gH::Int64`: 
 - `alLow::Float64`: 
 - `alTot::Float64`: 
 - `theta_f::Float64`: 
 - `theta_mid_neutral::Float64`: 
 - `al::Float64`: 
 - `be::Float64`: 
 - `B::Float64`: 
 - `bRange::Array{Float64,1}`:
 - `pposL::Float64`: 
 - `pposH::Float64`: 
 - `N::Int64`: 
 - `n::Int64`: 
 - `Lf::Int64`: 
 - `rho::Float64`: 
 - `TE::Float64`: 

"""
adap = parameters()

"""
	changeParameters()

Function to re-assign values to mutable struct *adap*. When values is not defined, it will be reset to the default value.

# Parameters
 - `gam_neg::Int64`: 
 - `gL::Int64`: 
 - `gH::Int64`: 
 - `alLow::Float64`: 
 - `alTot::Float64`: 
 - `theta_f::Float64`: 
 - `theta_mid_neutral::Float64`: 
 - `al::Float64`: 
 - `be::Float64`: 
 - `B::Float64`: 
 - `bRange::Array{Float64,1}`:
 - `pposL::Float64`: 
 - `pposH::Float64`: 
 - `N::Int64`: 
 - `n::Int64`: 
 - `Lf::Int64`: 
 - `rho::Float64`: 
 - `TE::Float64`: 

"""
function changeParameters(;gam_neg::Int64=-83,gL::Int64=10,gH::Int64=500,alLow::Float64=0.2,alTot::Float64=0.2,theta_f::Float64=1e-3,theta_mid_neutral::Float64=1e-3,al::Float64=0.184,be::Float64=0.000402,bRange::Array{Float64,1}=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.999],B::Float64=0.999,pposL::Float64=0.001,pposH::Float64=0.0,N::Int64=1000,n::Int64=661,Lf::Int64=10^6,rho::Float64=0.001,TE::Float64=5.0,diploid::Bool=true,convoluteBinomial::Bool=true)

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
	
	if diploid == false
		adap.NN = N
		adap.nn = n
	else
		adap.NN = (N*2)
		adap.nn = (n*2)
	end

	if convoluteBinomial == true
		adap.bn = Dict(bRange[i] => binomOp(bRange[i]) for i in 1:length(bRange))
	end
end

################################
###### Solving parameters ######
################################
"""
	Br(Lmax,theta)

Expected reduction in nucleotide diversity. Explored at [Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365):

```math
\\frac{\\pi}{\\pi_{0}} = e^{\\frac{4\\muL}{2rL+t}}
```

# Arguments
 - `Lmax::Int64`: non-coding flaking length
 - `theta::Float64`
# Returns
 - `Float64`: expected reduction in diversity given a non-coding length, mutation rate and defined recombination.
"""
function Br(Lmax::Int64,theta::Float64)

	ρ  	  = adap.rho
	t     = -1.0*adap.gam_neg/(adap.NN+0.0)
	μ     = theta/(2.0*adap.NN)
	r     = ρ/(2.0*adap.NN)

	return ℯ^(-4*μ*Lmax/(2*Lmax*r+t))
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value ) in a locus.
"""
	set_theta_f()

Find the optimum mutation given the expected reduction in nucleotide diversity (B value ) in a locus.

# Returns
 - `adap.theta_f::Float64`: changes adap.theta_f value.
"""
function set_theta_f()

	i(θ)         = Br(adap.Lf,θ)-adap.B
	theta_f      = Roots.find_zero(i,0.0)
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

"""
	setPpos()

Find the probabilty of positive selected alleles given the model. It solves a equation system taking into account fixations probabilities of weakly and strong beneficial alleles.

# Returns
 - `Tuple{Float64,Float64}`: weakly and strong beneficial alleles probabilites.
"""
function setPpos()
 	sc          = pyimport("scipy.optimize")
	pposL,pposH = sc.fsolve(solvEqns,(0.0,0.0))

	if pposL < 0.0
	 	pposL = 0.0
	end
	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
	if (pposH < 0.0 || pposH < 9e-15)
		pposH = 0.0
	end
	adap.pposL,adap.pposH = pposL, pposH
end

"""
	binomOp(B)

Site Frequency Spectrum convolution depeding on background selection values. Pass the SFS to a binomial distribution to sample the allele frequencies probabilites.

# Returns
 - `Array{Float64,2}`: convoluted SFS given a B value. It will be saved at *adap.bn*.
"""
function binomOp(B::Float64)

	NN2          = convert(Int64,ceil(adap.NN*adap.B))
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

"""

	phiReduction(gamma,ppos)


Reduction in fixation probabilty due to background selection and linkage. The formulas used have been subjected to several theoretical works ([Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365), [Hudson et al., 1995](https://www.genetics.org/content/141/4/1605), [Nordborg et al. 1995](https://doi.org/10.1017/S0016672300033619), [Barton NH., 1995](https://www.genetics.org/content/140/2/821)). 

The fixation probabilty of selected alleles are reduce by a factor ``\\phi``:

```math
\\phi(t,s) = e^{[\\frac{-2\\mu}{t(1+\\frac{rL}{t}+\\frac{2s}{t})}]}
```

Multiplying across all deleterious linkes sites, we find:

```math
\\Phi = \\prod_{1}^{L} = \\phi(t,s) = e^{\\frac{-2t\\mu(\\psi[1,\\frac{r+2s+L}{r}] - \\psi[1,\\frac{r(L+1)+2s+t}{r}])}{r^2}}
```

```math
\\phi(t,s) = e^{\\frac{-2t\\mu(\\psi[1,\\frac{r+2s+L}{r}] - \\psi[1,\\frac{r(L+1)+2s+t}{r}])}{r^2}}
```

# Arguments
 - `gamma::Int64`: selection coefficient.
	
# Returns
 - `Float64`: expected rate of positive fixations under background selection.

"""
function phiReduction(gammaValue::Int64)
	S  = abs(adap.gam_neg/(1.0*adap.NN))
	r  = adap.rho/(2.0*adap.NN)
	μ  = adap.theta_f/(2.0*adap.NN)
	s  = gammaValue/(adap.NN*1.0)

	Ψ0 = SpecialFunctions.polygamma(1,(s+S)/r)
	Ψ1 = SpecialFunctions.polygamma(1,(r+adap.Lf*r+s+S)/r)
	CC = 1.0

	return (ℯ^(-2.0*S*μ*(Ψ0-Ψ1)/(r^2)))
end

# function setPpos_nlsolve()

#    function f!(F,x)
#    	F[1] = alphaExpSimTot(x[1],x[2])-adap.alTot
#    	F[2] = alphaExpSimLow(x[1],x[2])-adap.alLow
#    end
   
#    pposL,pposH = nlsolve(f!,[0.0; 0.0]).zero

#    if pposL < 0.0
# 		pposL = 0.0
#    end

#    # Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
#    if (pposH < 0.0 || pposH < 9e-15)
# 	   pposH = 0.0
#    end

#    adap.pposL,adap.pposH = pposL, pposH
# end
