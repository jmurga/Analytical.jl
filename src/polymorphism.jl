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

Expected rate of neutral allele frequency reduce by background selection. The spectrum depends on the number of individual []

```math
\\mathbb{E}[Ps_{(x)}] = \\sum{x^{*}=x}{x^{*}=1}f_{B}(x)
```

# Return:
 - `Array{Float64}`: expected rate of neutral alleles frequencies.
"""
function DiscSFSNeutDown(param::parameters,binom::SparseMatrixCSC{Float64,Int64})

	NN2 = convert(Int64,ceil(param.NN*param.B))
	# Allocating variables

	neutralSfs(i::Int64) = 1.0/(i)

	x = collect(0:NN2)
	solvedNeutralSfs = neutralSfs.(x)
	replace!(solvedNeutralSfs, Inf => 0.0)

	# subsetDict = get(param.bn,param.B,1)
	# subsetDict = binom
	out::Array{Float64,1} = param.B*(param.thetaMidNeutral)*0.255*(binom*solvedNeutralSfs)
	# out = @view out[2:end-1]
	# out = out[2:end-1]

	return 	out
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
function DiscSFSSelPosDown(param::parameters,gammaValue::Int64,ppos::Float64,binom::SparseMatrixCSC{Float64,Int64})

	if ppos == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		exponentialType = Union{Float64,ArbFloat{48}}

		redPlus = phiReduction(param,gammaValue)

		# Solving sfs
		NN2 = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*param.B

		if isinf(exp(2*gammaCorrected))
			# Checked mpmath, BigFloat, DecFP.Dec128, Quadmath.Float128
			gammaExp1 = exp(ArbFloat(gammaCorrected*2,bits=24))
			gammaExp2 = exp(ArbFloat(gammaCorrected*-2,bits=24))
		else
			gammaExp1 = exp(gammaCorrected*2)
			gammaExp2 = exp(gammaCorrected*-2)
		end


		# Original
		# ppos*0.5*(ℯ^(2*gammaCorrected)*(1-ℯ^(-2.0*gammaCorrected*(1.0-i)))/((ℯ^(2*gammaCorrected)-1.0)*i*(1.0-i)))
		function positiveSfs(i::Float64,g1::T,g2::T,ppos::Float64) where {T<:Union{Float64,ArbFloat{48}}}
			if i > 0 && i < 1.0
				local tmp = ppos*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i)))
				return Float64(tmp)
			else
				return Float64(0.0)
			end
		end
		#
		#
		# function pSfs(x::Array{Float64,2},g1::T,g2::T,ppos::Float64) where {T<:Union{Float64,ArbFloat{48}}}
		# 	out = Array{Float64,1}[]
		# 	@inbounds @simd for i in x
		# 		if i > 0 && i < 1.0
		# 			tmp = ppos*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i)))
		# 		else
		# 			tmp = 0.0
		# 		end
		#
		# 	end
		# 	return out
		# end

		# Allocating outputs
		solvedPositiveSfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2,gammaExp1,gammaExp2,ppos))
		replace!(solvedPositiveSfs, NaN => 0.0)

		# subsetDict = get(param.bn,param.B,1)
		# out               = (param.thetaMidNeutral)*redPlus*0.745*(subsetDict*solvedPositiveSfs)
		out::Array{Float64,1} = (param.thetaMidNeutral)*redPlus*0.745*(binom*solvedPositiveSfs)
		# out = out[2:end-1]

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
# 		out               = (adap.thetaMidNeutral)*redPlus*0.745*(adap.bn[adap.B]*solvedPositiveSfs)
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
function DiscSFSSelNegDown(param::parameters,ppos::Float64,binom::SparseMatrixCSC{Float64,Int64})
	# subsetDict = get(param.bn,param.B,1)
	solvedNegative = DiscSFSSelNeg(param,ppos)
	# out::Array = param.B*(param.thetaMidNeutral)*0.745*(subsetDict*solvedNegative)
	out = param.B*(param.thetaMidNeutral)*0.745*(binom*solvedNegative)
	# out = @view out[2:end-1]

	# return out[2:end-1]
	return out

end

function DiscSFSSelNeg(param::parameters,ppos::Float64)

	beta     = param.be/(1.0*param.B)
	NN2      = convert(Int64, ceil(param.NN*param.B))
	xa       = collect(0:NN2)/NN2

	solveZ   = similar(xa)

	z(x::Float64,p::Float64=ppos) = (1.0-p)*(2.0^-param.al)*(beta^param.al)*(-SpecialFunctions.zeta(param.al,x+beta/2.0) + SpecialFunctions.zeta(param.al,(2+beta)/2.0))/((-1.0+x)*x)

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
	cumulativeSfs(sfsTemp)

Changing SFS considering all values above a frequency *x*. The original asymptotic-MK approach takes Pn(x) and Ps(x) as the number of polymorphic sites at frequency *x* rather than above *x*, but this approach scales poorly as sample size increases. We define the polymorphic spectrum as stated above since these quantities trivially have the same asymptote but are less affected by changing sample size.
"""
function cumulativeSfs(sfsTemp::Array)

	out      = Array{Float64}(undef, size(sfsTemp,1),size(sfsTemp,2))
	out[1,:] = sum(sfsTemp,dims=1)

	@simd for i in 2:(size(sfsTemp)[1])

		app = view(out,i-1,:) .- view(sfsTemp,i-1,:)

		if sum(app) > 0.0
			out[i,:] = app
		else
			out[i,:] = zeros(length(app))
		end
	end

	return out
end

"""
	reduceSfs(sfsTemp,bins)

Function to reduce the SFS into N bins.
"""
function reduceSfs(sfsTemp::Array,bins::Int64)

	freq  = collect(0:(size(sfsTemp,1)-1))/size(sfsTemp,1)
	h1    = fit(Histogram,freq,0:(1/bins):1)
	xmap1 = StatsBase.binindex.(Ref(h1), freq)

	tmp = hcat(sfsTemp,xmap1)
	out = Array{Float64}(undef,bins-1	,size(sfsTemp,2))
	vIter =  convert(Array,unique(xmap1)')
	@simd for i = eachindex(vIter)
		@inbounds out[i,:] = sum(tmp[tmp[:,end] .== i,1:end-1],dims=1)
	end

	return (out')
end
