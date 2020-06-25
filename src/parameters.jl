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
	diploid::Bool			   = false

	NN::Int64 = 2*N
	nn::Int64 = 2*n

	bn::Dict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()

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
# global adap = parameters()

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
function Br(param::parameters,theta::Float64)

	ρ::Float64     = param.rho
	t::Float64     = -1.0*param.gam_neg/(param.NN+0.0)
	μ::Float64     = theta/(2.0*param.NN)
	r::Float64     = ρ/(2.0*param.NN)

	out::Float64  = ℯ^(-4*μ*param.Lf/(2*param.Lf*r+t))
	return out
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value ) in a locus.
"""
	set_theta_f(param)

Find the optimum mutation given the expected reduction in nucleotide diversity (B value) in a locus.

# Returns
 - `adap.theta_f::Float64`: changes adap.theta_f value.
"""
function set_theta_f(param::parameters)

	i(θ,p=param) = Br(p,θ)-p.B
	theta_f      = Roots.find_zero(i,0.0)
	param.theta_f = theta_f
end

function alphaExpSimLow(param::parameters,pposL::Float64,pposH::Float64)
	return fixPosSim(param,param.gL,0.5*pposL)/(fixPosSim(param,param.gL,0.5*pposL)+fixPosSim(param,param.gH,0.5*pposH) + fixNegB(param,0.5*pposL+0.5*pposH))
end

function alphaExpSimTot(param::parameters,pposL::Float64,pposH::Float64)
	return (fixPosSim(param,param.gL,0.5*pposL)+fixPosSim(param,param.gH,0.5*pposH))/(fixPosSim(param,param.gL,0.5*pposL)+fixPosSim(param,param.gH,0.5*pposH)+fixNegB(param,0.5*pposL+0.5*pposH))
end

function solvEqns(param,config)

	pposL,pposH = config
	return (alphaExpSimTot(param,pposL,pposH)-param.alTot,alphaExpSimLow(param,pposL,pposH)-param.alLow)
end

"""
	setPpos(param)

Find the probabilty of positive selected alleles given the model. It solves a equation system taking into account fixations probabilities of weakly and strong beneficial alleles.

# Returns
 - `Tuple{Float64,Float64}`: weakly and strong beneficial alleles probabilites.
"""

function setPpos(param::parameters)

	function f!(F,x,param=param)
		F[1] = alphaExpSimTot(param,x[1],x[2])-param.alTot
		F[2] = alphaExpSimLow(param,x[1],x[2])-param.alLow
	end
 
	pposL,pposH = NLsolve.nlsolve(f!,[0.0; 0.0]).zero
 
	if pposL < 0.0
		 pposL = 0.0
	end  
	if pposH < 0.0
		 pposH = 0.0
	end

	param.pposL,param.pposH = pposL, pposH

 end

 
# function setPpos()
#  	sc          = pyimport("scipy.optimize")
# 	pposL,pposH = sc.fsolve(solvEqns,(0.0,0.0))

# 	if pposL < 0.0
# 	 	pposL = 0.0
# 	end
# 	if pposH < 0.0
# 		pposH = 0.0
#    end
# 	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results

# 	param.pposL,param.pposH = pposL, pposH
# end

"""
	binomOp(param)

Site Frequency Spectrum convolution depeding on background selection values. Pass the SFS to a binomial distribution to sample the allele frequencies probabilites.

# Returns
 - `Array{Float64,2}`: convoluted SFS given a B value. It will be saved at *adap.bn*.
"""
function binomOp(param::parameters)
        
    bn = Dict(param.bRange[i] => zeros(param.nn+1,param.NN) for i in 1:length(param.bRange))

    for bVal in param.bRange

        NN2          = convert(Int64,ceil(param.NN*bVal))
        samples      = collect(0:param.nn)
        samplesFreqs = collect(0:NN2)
        samplesFreqs = permutedims(samplesFreqs/NN2)
    
        f(x) = Distributions.Binomial(param.nn,x)
        z    = f.(samplesFreqs)     
    
        out  = Distributions.pdf.(z,samples)
		outS  = round.(out,digits=10)
        param.bn[bVal] = SparseArrays.sparse(outS)

	end
end


"""

	phiReduction(param,gamma,ppos)


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
function phiReduction(param::parameters,gammaValue::Int64)
	S::Float64 = abs(param.gam_neg/(1.0*param.NN))
	r::Float64 = param.rho/(2.0*param.NN)
	μ::Float64 = param.theta_f/(2.0*param.NN)
	s::Float64 = gammaValue/(param.NN*1.0)

	Ψ0::Float64 = SpecialFunctions.polygamma(1,(s+S)/r)
	Ψ1::Float64 = SpecialFunctions.polygamma(1,(r+param.Lf*r+s+S)/r)
	CC::Float64 = 1.0

	out::Float64 = (ℯ^(-2.0*S*μ*(Ψ0-Ψ1)/(r^2)))
	return out
end
