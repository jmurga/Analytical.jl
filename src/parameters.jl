################################
####### Define parameters ######
################################
import Parameters: @with_kw

"""
	parameters(
		gamNeg::Int64,
		gL::Int64,
		gH::Int64,
		alLow::Float64,
		alTot::Float64,
		thetaF::Float64,
		thetaMidNeutral::Float64,
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
 - `gamNeg::Int64`:
 - `gL::Int64`:
 - `gH::Int64`:
 - `alLow::Float64`:
 - `alTot::Float64`:
 - `thetaF::Float64`:
 - `thetaMidNeutral::Float64`:
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
@with_kw mutable struct parameters
	gamNeg::Int64              = -457
	gL::Int64                  = 10
	gH::Int64                  = 500
	alLow::Float64             = 0.2
	alTot::Float64             = 0.4
	thetaF::Float64            = 1e-3
	thetaMidNeutral::Float64   = 1e-3
	al::Float64                = 0.184
	be::Float64                = 0.000402
	B::Float64                 = 0.999 
	bRange::Array{Float64,2}   = permutedims(push!(collect(0.1:0.025:0.975),0.999))
	pposL::Float64             = 0
	pposH::Float64             = 0.001
	N::Int64                   = 1000
	n::Int64                   = 500
	Lf::Int64                  = 2*10^5
	rho::Float64               = 0.001
	TE::Float64                = 5.0
	diploid::Bool              = false

	NN::Int64 = 2*N
	nn::Int64 = 2*n
	dac::Array{Int64,1} = [2,4,5,10,20,50,200,500,700]

	#=bn::Dict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()=#
	#=neut::Dict = Dict{Float64,Array{Float64,1}}()=#
end

@with_kw mutable struct binomialDict
	bn::Dict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()
end
#=const binomialDict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()=#

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
	t::Float64     = -1.0*param.gamNeg/(param.NN+0.0)
	μ::Float64     = theta/(2.0*param.NN)
	r::Float64     = ρ/(2.0*param.NN)

	out::Float64  = ℯ^(-4*μ*param.Lf/(2*param.Lf*r+t))
	return out
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value ) in a locus.
"""
	setThetaF!(param)

Find the optimum mutation given the expected reduction in nucleotide diversity (B value) in a locus.

# Returns
 - `adap.thetaF::Float64`: changes adap.thetaF value.
"""
function setThetaF!(param::parameters)

	i(θ,p=param) = Br(p,θ)-p.B
	thetaF      = Roots.find_zero(i,0.0)
	param.thetaF = thetaF
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
	setPpos!(param)

Find the probabilty of positive selected alleles given the model. It solves a equation system taking into account fixations probabilities of weakly and strong beneficial alleles.

# Returns
 - `Tuple{Float64,Float64}`: weakly and strong beneficial alleles probabilites.
"""

function setPpos!(param::parameters)

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

"""
	binomOp(param)

Site Frequency Spectrum convolution depeding on background selection values. Pass the SFS to a binomial distribution to sample the allele frequencies probabilites.

# Returns
 - `Array{Float64,2}`: convoluted SFS given a B value. It will be saved at *adap.bn*.
"""
function binomOp!(param::parameters,convolutedBn::Dict)

	#=bn = Dict(param.bRange[i] => spzeros(param.nn+1,param.NN) for i in 1:length(param.bRange))=#

	for bVal in param.bRange

		NN2          = convert(Int64,ceil(param.NN*bVal))
		samples      = collect(1:(param.nn-1))
		pSize        = collect(0:NN2)
		samplesFreqs = permutedims(pSize/NN2)
		neutralSfs   = @. 1/pSize
		replace!(neutralSfs, Inf => 0.0)


		f(x) = Distributions.Binomial(param.nn,x)
		z    = f.(samplesFreqs)

		out  = Distributions.pdf.(z,samples)
		out  = round.(out,digits=10)
		outS = SparseArrays.dropzeros(SparseArrays.sparse(out))
		convolutedBn[bVal] = outS
		# param.bn[bVal] = outS
		#=param.neut[bVal] = round.(param.B*(param.thetaMidNeutral)*0.25*(outS*neutralSfs),digits=10)=#

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
	S::Float64 = abs(param.gamNeg/(1.0*param.NN))
	r::Float64 = param.rho/(2.0*param.NN)
	μ::Float64 = param.thetaF/(2.0*param.NN)
	s::Float64 = gammaValue/(param.NN*1.0)

	Ψ0::Float64 = SpecialFunctions.polygamma(1,(s+S)/r)
	Ψ1::Float64 = SpecialFunctions.polygamma(1,(r+param.Lf*r+s+S)/r)
	CC::Float64 = 1.0

	out::Float64 = (ℯ^(-2.0*S*μ*(Ψ0-Ψ1)/(r^2)))
	return out
end

"""
	alphaByFrequencies(gammaL,gammaH,pposL,param.pposH,data,nopos)

Analytical α(x) estimation. Solve α(x) from the expectation generally. We used the expected rates of divergence and polymorphism to approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region.

```math
\\mathbb{E}[\\alpha_{x}] =  1 - \\left(\\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{N}]}\\frac{\\mathbb{E}[P_{N}]}{\\mathbb{E}[P_{S}]}\\right)
```

# Arguments
 - `gammaL::Int64`: strength of weakly positive selection
 - `gammaH::Int64`: strength of strong positive selection
 - `pposL`::Float64: probability of weakly selected allele
 - `param.pposH`::Float64: probability of strong selected allele
 - `data::Array{Any,1}`: Array containing the total observed divergence, polymorphism and site frequency spectrum.
 - `nopos::String("pos","nopos","both")`: string to perform α(x) account or not for both positive selective alleles.

# Returns
 - `Array{Float64,1}` α(x).
"""
function analyticalAlpha(;param::parameters,convolutedSamples::binomialDict)

	##############################################################
						# Solve the model  #
	##############################################################
	B = param.B

	setThetaF!(param)
	thetaF = param.thetaF
	# Solve the probabilities of fixations without background selection
	## First set non-bgs
	param.B = 0.999
	## Solve the mutation rate
	setThetaF!(param)
	## Solve the probabilities
	setPpos!(param)
	# Return to the original values
	param.thetaF = thetaF
	param.B = B

	##############################################################
	# Accounting for positive alleles segregating due to linkage #
	##############################################################

	# Fixation
	fN     = param.B*fixNeut(param)
	fNeg   = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
	fPosL  = fixPosSim(param,param.gL,0.5*param.pposL)
	fPosH  = fixPosSim(param,param.gH,0.5*param.pposH)

	ds = fN
	dn = fNeg + fPosL + fPosH

	## Polymorphism
	neut = DiscSFSNeutDown(param,convolutedSamples.bn[param.B])

	selH = DiscSFSSelPosDown(param,param.gH,param.pposH,convolutedSamples.bn[param.B])
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL,convolutedSamples.bn[param.B])
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,convolutedSamples.bn[param.B])
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN),false)

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	# ps = @. neut / (sel+neut)
	# pn = @. sel / (sel+neut)

	## Outputs
	α = @. 1 - ((ds/dn) * (sel/neut))


	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	fN_nopos     = fN*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fNeg_nopos   = fNeg*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fPosL_nopos  = fPosL*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fPosH_nopos  = fPosH*(param.thetaMidNeutral/2.)*param.TE*param.NN

	ds_nopos = fN_nopos
	dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

	## Polymorphism
	sel_nopos = selN
	ps_nopos = @. neut / (sel_nopos + neut)
	pn_nopos = @. sel_nopos / (sel_nopos + neut)

	α_nopos = 1 .- (ds_nopos/dn_nopos) .* (sel_nopos./neut)

	##########
	# Output #
	##########
	return (α,α_nopos)
end
