################################
####### Define parameters ######
################################
"""
Mutable structure containing the variables required to solve the analytical approach. All the functions are solve using the internal values of the structure. You should declare a mutable structure to the perform the analytical estimations.

# Parameters
 - `gamNeg::Int64`: Selection coefficient for deleterious alleles
 - `gL::Int64`: Selection coefficient for weakly benefical alleles
 - `gH::Int64`: Selection coefficient for strongly benefical alleles
 - `alLow::Float64`: Proportion of α due to weak selection
 - `alTot::Float64`: α
 - `thetaF::Float64`: Mutation rate defining BGS strength
 - `thetaMidNeutral::Float64`: Mutation rate on coding region
 - `al::Float64`: DFE shape parameter 
 - `be::Float64`: DFE scale parameter
 - `B::Float64`: BGS strength
 - `B_bins::Array{Float64,1}`: BGS values to simulate
 - `pposL::Float64`: Fixation probabily of weakly beneficial alleles
 - `pposH::Float64`: Fixation probabily of strongly beneficial alleles
 - `N::Int64`: Population size
 - `n::Int64`: Sample size
 - `Lf::Int64`: Flanking region length
 - `rho::Float64`: Recombination rate
 - `TE::Float64`

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
	B_bins::Array{Float64,1}   = push!(collect(0.1:0.025:0.975),0.999)
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
	θᵣ::Array{Float64,1} = fill(thetaMidNeutral,nn-1)
	dac::Array{Int64,1} = [2,4,5,10,20,50,200,500,700]

	binom::Dict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()
end

"""
Mutable structure containing the downsampled SFS. 

# Returns
 - `bn::Dict`: SparseMatrixCSC containing the binomial convolution

"""
#=@with_kw mutable struct binomial_dict
	bn::Dict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()
end=#
#=const binomial_dict = Dict{Float64,SparseMatrixCSC{Float64,Int64}}()=#

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
 - `param::parameters`
 - `theta::Float64`
# Returns
 - `Float64`: expected reduction in diversity given a non-coding length, mutation rate and defined recombination.
"""
function Br(param::parameters,theta::Float64)

	@unpack rho, gamNeg, NN, NN, NN, Lf, Lf = param;

	ρ::Float64     = rho
	t::Float64     = -1.0*gamNeg/(NN+0.0)
	μ::Float64     = theta/(2.0*NN)
	r::Float64     = ρ/(2.0*NN)

	out::Float64   = ℯ^(-4*μ*Lf/(2*Lf*r+t))
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

	@unpack B = param
	i(θ,p=param,b=B) = Br(p,θ) - B;
	thetaF           = find_zero(i,0.0);
	@pack! param    = thetaF;
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

	pposL,pposH = nlsolve(f!,[0.0; 0.0]).zero

	if pposL < 0.0
		 pposL = 0.0
	end
	
	if pposH < 0.0
		 pposH = 0.0
	end

	@pack! param = pposL, pposH;

end

"""
	binomOp(param)

Binomial convolution to sample the allele frequencies probabilites depending on background selection values and sample size.

# Arguments
 - `param::parameters`
 - `convoluted_samples::binomial_dict`

# Returns
 - `Array{Float64,2}`: convoluted SFS for each B value defined in the model (param.B_bins). The estimations are saved at *convolutedBn.bn*.
"""
function binomOp!(param::parameters)

	# bn = Dict(param.B_bins[i] => spzeros(param.nn+1,param.NN) for i in 1:length(param.B_bins))

	for b_val in param.B_bins

		NN2          = convert(Int64,ceil(param.NN*b_val))
		samples      = collect(1:(param.nn-1))
		pSize        = collect(0:NN2)
		samplesFreqs = permutedims(pSize/NN2)
		neutralSfs   = @. 1/pSize
		replace!(neutralSfs, Inf => 0.0)


		f(x,y=param.nn) = Binomial(y,x)
		z    = f.(samplesFreqs)

		out  = pdf.(z,samples)
		out  = round.(out,digits=10)
		out_S = SparseArrays.dropzeros(SparseArrays.sparse(out))
		param.binom[b_val] = out_S
		# bn[b_val] = out_S
		#=param.neut[bVal] = round.(param.B*(param.thetaMidNeutral)*0.25*(outS*neutralSfs),digits=10)=#

	end
end

function binomOp!(NN::Int64,nn::Int64,B_bins::Vector{Float64})

	bn = Dict(B_bins[i] => spzeros(nn+1,NN) for i in 1:length(B_bins))

	for b_val in B_bins

		NN2          = convert(Int64,ceil(NN*b_val))
		samples      = collect(1:(nn-1))
		pSize        = collect(0:NN2)
		samplesFreqs = permutedims(pSize/NN2)
		neutralSfs   = @. 1/pSize
		replace!(neutralSfs, Inf => 0.0)


		f(x,y=param.nn) = Binomial(y,x)
		z    = f.(samplesFreqs)

		out  = pdf.(z,samples)
		out  = round.(out,digits=10)
		out_S = SparseArrays.dropzeros(SparseArrays.sparse(out))
		
		bn[b_val] = out_S

	end
	return(bn)
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

	Ψ0::Float64 = polygamma(1,(s+S)/r)
	Ψ1::Float64 = polygamma(1,(r+param.Lf*r+s+S)/r)
	CC::Float64 = 1.0

	out::Float64 = (ℯ^(-2.0*S*μ*(Ψ0-Ψ1)/(r^2)))
	return out
end

"""
	analytical_alpha(param, convoluted_samples)

Analytical α(x) estimation. Solve α(x) generally. We used the expected rates of divergence and polymorphism to approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region.

```math
\\mathbb{E}[\\alpha_{x}] =  1 - \\left(\\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{N}]}\\frac{\\mathbb{E}[P_{N}]}{\\mathbb{E}[P_{S}]}\\right)
```

# Arguments
 - `param::parameters`
 - `convoluted_samples::binomial_dict`
# Returns
 - `Array{Float64,1}` α(x).
"""
function analytical_alpha(;param::parameters)

	################################################################
		# Solve the model similarly to original python mktest  #	
	################################################################
	B = param.B

	setThetaF!(param)
	thetaF = param.thetaF
	# Solve the probabilities of fixations without BGS
	## Set non-bgs
	param.B = 0.999
	## Solve the mutation rate
	setThetaF!(param)
	## Solve the probabilities
	setPpos!(param)
	# Return to the original values
	param.thetaF = thetaF
	param.B = B

	################################################################
	 # Accounting for positive alleles segregating due to linkage # 
	################################################################

	# Fixation
	fN     = param.B*fixNeut(param)
	fNeg   = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
	fPosL  = fixPosSim(param,param.gL,0.5*param.pposL)
	fPosH  = fixPosSim(param,param.gH,0.5*param.pposH)

	ds = fN
	dn = fNeg + fPosL + fPosH

	## Polymorphism
	neut = DiscSFSNeutDown(param)

	selH = if isinf(exp(param.gH * 2))
		DiscSFSSelPosDownArb(param,param.gH,param.pposH,)
	else
		DiscSFSSelPosDown(param,param.gH,param.pposH)
	end

	selL = DiscSFSSelPosDown(param,param.gL,param.pposL)	
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL)

	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	tmp = cumulative_sfs(hcat(neut,selH,selL,selN),false)

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	# ps = @. neut / (sel+neut)
	# pn = @. sel / (sel+neut)

	## Outputs
	α = @. 1 - ((ds/dn) * (sel/neut))


	#############################################################
	#		Accounting for neutral and deleterious alleles		#	#############################################################
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
