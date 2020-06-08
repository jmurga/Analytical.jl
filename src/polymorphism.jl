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
function DiscSFSNeutDown()

	NN2 = convert(Int64,ceil(adap.NN*adap.B))
	# Allocating variables
	x                = Array{Int64}(undef,NN2 + 1)
	solvedNeutralSfs = Array{Float64}(undef,NN2 + 1)
	out              = Array{Float64}(undef,NN2 + 1)

	function neutralSfs(i::Int64)
		if i > 0 && i < NN2
			 return 1.0/(i)
		end
		return 0.0
	end

	x                = collect(0:NN2)
	solvedNeutralSfs = x .|> neutralSfs

	# c = similar(adap.bn[adap.B][:,1])

	out::Array{Float64} = adap.B*(adap.theta_mid_neutral)*0.255*(adap.bn[adap.B]*solvedNeutralSfs)

	return 	view(out,2:lastindex(out)-1,:)
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
function DiscSFSSelPosDown(gammaValue::Int64,ppos::Float64)

	if ppos == 0.0
		out = zeros(Float64,adap.nn + 1)
	else

		red_plus = phiReduction(gammaValue)

		# Solving sfs
		NN2 = convert(Int64,ceil(adap.NN*adap.B))
		xa  = collect(0:NN2)
		xa  = xa/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*adap.B
		if isinf(exp(2*gammaCorrected))
			# Checked mpmath, BigFloat, DecFP.Dec128, Quadmath.Float128
			gammaExp1 = exp(ArbFloat(gammaCorrected*2))
			gammaExp2 = exp(ArbFloat(gammaCorrected*-2))
		else
			gammaExp1 = exp(gammaCorrected*2)
			gammaExp2 = exp(gammaCorrected*-2)
		end

		# Original
		# ppos*0.5*(ℯ^(2*gammaCorrected)*(1-ℯ^(-2.0*gammaCorrected*(1.0-i)))/((ℯ^(2*gammaCorrected)-1.0)*i*(1.0-i)))
		function positiveSfs(i,gammaExp1=gammaExp1,gammaExp2=gammaExp2,ppos=ppos)
			if i > 0 && i < 1.0
				# if gammaExp1 is Float64
				# if gammaValue < 20
				# 	randGammaExp = rand(Distributions.Exponential(gammaCorrected))
				# 	gammaExp1    = exp(randGammaExp*2)
				# 	gammaExp2    = exp(randGammaExp*-2)
				# end
				local out = ppos*0.5*(gammaExp1*(1- gammaExp2^(1.0-i))/((gammaExp1-1.0)*i*(1.0-i)))
				return Float64(out)
			else
				return 0.0
			end
		end

		# Allocating outputs
		solvedNeutralSfs = Array{Float64}(undef,NN2 + 1)
		out              = Array{Float64}(undef,NN2 + 1)

		solvedPositiveSfs = (1.0/(NN2)) * (positiveSfs.(xa))
		replace!(solvedPositiveSfs, NaN => 0.0)
		out               = (adap.theta_mid_neutral)*red_plus*0.745*(adap.bn[adap.B]*solvedPositiveSfs)
	end

	return view(out,2:lastindex(out)-1,:)
end

# function DiscSFSSelPosDown(gammaValue::Int64,ppos::Float64)

# 	if ppos == 0.0
# 		out = zeros(Float64,adap.nn + 1)
# 	else

# 		red_plus = phiReduction(gammaValue)

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
# 		out               = (adap.theta_mid_neutral)*red_plus*0.745*(adap.bn[adap.B]*solvedPositiveSfs)
# 	end

# 	return view(out,2:lastindex(out)-1,:)
# end

######Slightly deleterious######
"""

	DiscSFSSelNegDown(ppos)

Expected rate of positive selected allele frequency reduce by background selection. Spectrum drawn on a gamma DFE. It depends on the number of individuals.

# Arguments
 - `ppos::Float64`: positive selected alleles probabilty.
# Return:
 - `Array{Float64}`: expected negative selected alleles frequencies.
"""
function DiscSFSSelNegDown(ppos::Float64)
	out::Array = adap.B*(adap.theta_mid_neutral)*0.745*(adap.bn[adap.B]*DiscSFSSelNeg(ppos))
	return out[2:lastindex(out)-1]
end

function DiscSFSSelNeg(ppos::Float64)

	beta     = adap.be/(1.0*adap.B)
	NN2      = convert(Int64, ceil(adap.NN*adap.B))
	xa       = collect(0:NN2)/NN2

	solveZ   = similar(xa)

	z(x::Float64,ppos::Float64=ppos) = (1.0-ppos)*(2.0^-adap.al)*(beta^adap.al)*(-SpecialFunctions.zeta(adap.al,x+beta/2.0) + SpecialFunctions.zeta(adap.al,(2+beta)/2.0))/((-1.0+x)*x)

	solveZ   = xa .|> z

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
function cumulativeSfs(sfsTemp)

	out      = Array{Float64}(undef, size(sfsTemp,1) + 1,size(sfsTemp,2))
	out[1,:] = sum(sfsTemp,dims=1)

	for i in 2:(size(sfsTemp)[1]+1)

		app = out[i-1,:] .- sfsTemp[i-1,:]

		if sum(app) > 0.0
			out[i,:] = app
		else
			out[i,:] = zeros(length(app))
		end
	end

	return out
end
