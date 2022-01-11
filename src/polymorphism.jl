################################
######    Polymorphism    ######
################################
# Expected number of polymorphism above frequency x from the standard diffusion theory
# f(x) = ∫(s) θs* (1/x(1-x)) * ( (ℯ^(Ns)*(1-ℯ^(-4Ns(1-x))) ) / (ℯ^(4Ns-1)))
# Convolote with the binomial to obtain the downsampled sfs
# E[P(x)] = ∑_x=x*→x=1 fB(x*)
############Neutral#############

"""

	DiscSFSNeutDown()

Expected rate of neutral allele frequency reduce by backgrou	nd selection. The spectrum depends on the number of individual []

```math
\\mathbb{E}[Ps_{(x)}] = \\sum{x^{*}=x}{x^{*}=1}f_{B}(x)
```

# Return:
 - `Array{Float64}`: expected rate of neutral alleles frequencies.
"""
function DiscSFSNeutDown(param::parameters)

	NN2::Int64 = convert(Int64,ceil(param.NN*param.B))
	# Allocating variables

	neutralSfs(i::Int64) = 1.0/(i)

	x::Vector{Int64} = collect(0:NN2)
	solvedNeutralSfs::Vector{Float64} = neutralSfs.(x)
	replace!(solvedNeutralSfs, Inf => 0.0)

	# subsetDict = get(param.binom,param.B,1)
	# subsetDict = binom
	out::Array{Float64,1} = param.B * (param.thetaMidNeutral) * 0.25 * (param.binom[param.B]*solvedNeutralSfs)
	#=out::Array{Float64,1} = param.B .* (param.θᵣ) .*0.25 .* (binom*solvedNeutralSfs)=#
	return out
end

############Positive############
# Variable gamma in function changed to gammaValue to avoid problem with exported SpecialFunctions.gamma
"""

	DiscSFSSelPosDown(gammaValue,ppos)

Expected rate of positive selected allele frequency reduce by background selection. The spectrum depends on the number of individuals.

# Arguments
 - `gammaValue::Int64`: selection strength.
 - `ppos::Float64`: positive selected alleles probabilty.
# Return:
 - `Array{Float64}`: expected positive selected alleles frequencies.
"""
function DiscSFSSelPosDown(param::parameters,gammaValue::Int64,ppos::Float64)

	if ppos == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		redPlus = phiReduction(param,gammaValue)

		# Solving sfs
		NN2 = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*param.B

		gammaExp1 = exp(gammaCorrected*2)
		gammaExp2 = exp(gammaCorrected*-2)

		positiveSfs(i::Float64,g1::Float64=gammaExp1,g2::Float64=gammaExp2,ppos::Float64=ppos) = Float64(ppos*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i))))

		# Original
		# ppos*0.5*(ℯ^(2*gammaCorrected)*(1-ℯ^(-2.0*gammaCorrected*(1.0-i)))/((ℯ^(2*gammaCorrected)-1.0)*i*(1.0-i)))

		# Allocating outputs
		solvedPositiveSfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2))
		replace!(solvedPositiveSfs, NaN => 0.0)

		# subsetDict = get(param.binom,param.B,1)
		# out               = (param.thetaMidNeutral)*redPlus*0.75*(subsetDict*solvedPositiveSfs)
		out::Array{Float64,1} = param.thetaMidNeutral * redPlus * 0.75 * (param.binom[param.B]*solvedPositiveSfs)

	end

	return out
end

function DiscSFSSelPosDownArb(param::parameters,gammaValue::Int64,ppos::Float64)

	if ppos == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		redPlus = phiReduction(param,gammaValue)

		# Solving sfs
		NN2  = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*param.B
		gammaExp1::Quadmath.Float128 = exp(
			Quadmath.Float128(gammaCorrected*2))
		gammaExp2::Quadmath.Float128 = exp(Quadmath.Float128(gammaCorrected*-2))

		positiveSfs(i::Float64,g1::Quadmath.Float128=gammaExp1,g2::Quadmath.Float128=gammaExp2,ppos::Float64=ppos) = Float64(ppos*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i))))
		# Allocating outputs
		solvedPositiveSfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2))
		replace!(solvedPositiveSfs, NaN => 0.0)
		out::Array{Float64,1} = param.thetaMidNeutral * redPlus * 0.75 * (param.binom[param.B]*solvedPositiveSfs)
		#=out::Array{Float64,1} = param.θᵣ .* redPlus .* 0.75 .* (binom*solvedPositiveSfs)=#

	end

	return out
end

# function DiscSFSSelPosDown(gammaValue::Int64,ppos::Float64)

# 	if ppos == 0.0
# 		out = zeros(Float64,adap.nn + 1)
# 	else

# 		redPlus = phiReduction(gammaValue)

# 		# Solving sfs
# 		NN2 = convert(Int64,ceil(adap.NN*adap.B))
# 		xa  = collect(0:NN2)
# 		xa  = xa/(NN2)

		# function positiveSfs(i,gammaCorrected=gammaValue*adap.B,ppos=ppos)
		# 	if i > 0 && i < 1.0
		# 		return ppos*0.5*(
		# 			ℯ^(2*gammaCorrected)*(1-ℯ^(-2.0*gammaCorrected*(1.0-i)))/((ℯ^(2*gammaCorrected)-1.0)*i*(1.0-i)))
		# 	end
		# 	return 0.0
		# end

# 		# Allocating outputs
# 		solvedNeutralSfs = Array{Float64}(undef,NN2 + 1)
# 		out              = Array{Float64}(undef,NN2 + 1)

# 		solvedPositiveSfs = (1.0/(NN2)) * (xa .|> positiveSfs)
# 		replace!(solvedPositiveSfs, NaN => 0.0)
# 		out               = (adap.thetaMidNeutral)*redPlus*0.75*(adap.bn[adap.B]*solvedPositiveSfs)
# 	end

# 	return view(out,2:lastindex(out)-1,:)
# end

######Slightly deleterious######
"""

	DiscSFSSelNegDown(param,ppos)

Expected rate of positive selected allele frequency reduce by background selection. Spectrum drawn on a gamma DFE. It depends on the number of individuals.

# Arguments
 - `ppos::Float64`: positive selected alleles probabilty.
# Return:
 - `Array{Float64}`: expected negative selected alleles frequencies.
"""
function DiscSFSSelNegDown(param::parameters,ppos::Float64)
	# subsetDict = get(param.binom,param.B,1)
	solvedNegative = DiscSFSSelNeg(param,ppos)
	out = param.B .* (param.thetaMidNeutral) .* 0.75 .* (param.binom[param.B]*solvedNegative)
	#=out = param.B .* (param.θᵣ) .* 0.75 .* (binom*solvedNegative)=#

	return out
end

function DiscSFSSelNeg(param::parameters,ppos::Float64)

	beta     = param.be/(1.0*param.B)
	NN2      = convert(Int64, ceil(param.NN*param.B))
	xa       = collect(0:NN2)/NN2

	solveZ   = similar(xa)

	z(x::Float64,p::Float64=ppos) = (1.0-p)*(2.0^-param.al)*(beta^param.al)*(-zeta(param.al,x+beta/2.0) + zeta(param.al,(2+beta)/2.0))/((-1.0+x)*x)

	solveZ   = z.(xa)

	if (solveZ[1] == Inf || isnan(solveZ[1]))
		solveZ[1] = 0.0
	end
	if (solveZ[lastindex(solveZ)] == Inf || isnan(solveZ[lastindex(solveZ)]))
		solveZ[lastindex(solveZ)] = 0.0
	end

	return 1.0/(NN2+0.0).*solveZ
end

"""
	cumulative_sfs(sfs_tmp)

Changing SFS considering all values above a frequency *x*. The original asymptotic-MK approach takes Pn(x) and Ps(x) as the number of polymorphic sites at frequency *x* rather than above *x*, but this approach scales poorly as sample size increases. We define the polymorphic spectrum as stated above since these quantities trivially have the same asymptote but are less affected by changing sample size.
"""
function cumulative_sfs(sfs_tmp::Array,freqs::Bool=true)

	out      = Array{Float64}(undef, size(sfs_tmp,1),size(sfs_tmp,2))

	if freqs
		idx = 2
	else
		idx = 1
	end

	out[1,idx:end] = sum(sfs_tmp[:,idx:end],dims=1)

	@simd for i in 2:(size(sfs_tmp)[1])

		#=app = view(out,i-1,:) .- view(sfs_tmp,i-1,:)=#
		app = out[i-1,idx:end] .- sfs_tmp[i-1,idx:end]

		if sum(app) > 0.0
			out[i,idx:end] = app
		else
			out[i,idx:end] = zeros(length(app))
		end
	end

	if freqs
		out[:,1] = sfs_tmp[:,1]
	end
	
	return out
end

"""
	reduce_sfs(sfs_tmp,bins)

Function to bin the SFS into a sample of N individuals.
"""
function reduce_sfs(sfs_tmp::Array,bins::Int64)

	n   = (bins*2) - 1 
	f   = sfs_tmp[:,1]
	sfs = sfs_tmp[:,2:end]
	
	b    = collect(1/n:1/n:1)
	inds = searchsortedfirst.(Ref(b), f)
	out  = zeros((n,size(sfs_tmp,2)))
	out[:,1] = unique(inds)
	sfs_grouped = hcat(inds,sfs)
	
	for i in unique(inds)
		out[out[:,1] .== i,2:end] = sum(sfs_grouped[sfs_grouped[:,1] .== i,2:end],dims=1)
	end
	out[:,1] = b

	return(out)
end
