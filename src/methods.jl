"""
	asympFit(x,columns)

Function to estimate the asmpytotic value of α(x).

# Arguments
 - `alphaValues::Array{Float64,1}`: α(x) array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function asympFit(alphaValues::Array{Float64,1})

	# Model
	asympModel(x,p) = @. p[1] + p[2]*exp(-x*p[3])

	# Fit values
	fitted1    = LsqFit.curve_fit(asympModel,collect(1:size(alphaValues,1)),alphaValues,[-1.0,-1.0,1.0];lower=[-1.0,-1.0,1.0],upper=[1.0, 1.0, 10.0])
	fitted2    = LsqFit.curve_fit(asympModel,collect(1:size(alphaValues,1)),alphaValues,fitted1.param)
	asymp      = asympModel(size(alphaValues,1),fitted2.param)

	ciLow, ciHigh   = try
		LsqFit.confidence_interval(fitted2)[1][1],LsqFit.confidence_interval(fitted2)[1][2]
	catch err
		(0.0,0.0)
	end

	return asymp,[ciLow,ciHigh],fitted2.param
	# return asymp
end

"""
	imputedMK(sfs,divergence,m,cutoff)

Function to estimate the asmpytotic value of α(x).

# Arguments
 - `alphaValues::Array{Float64,1}`: α(x) array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function impMK(;sfs::Array,divergence::Array,m::T,cutoff::Float64=0.15) where {T<:Union{Nothing,Array}}

	output = OrderedDict{String,Float64}()

	ps = sum(sfs[:,2])
	pn = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	
	deleterious = 0
	### Estimating slightly deleterious with pn/ps ratio
	fltLow = (sfs[:, 1] .<= cutoff)
	pnLow   = sum(sfs[fltLow,2])
	psLow   = sum(sfs[fltLow,3])

	fltInter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
	pnInter = sum(sfs[fltInter,2])
	psInter = sum(sfs[fltInter,3])

	ratioPs       = psLow / psInter
	deleterious   = pnLow - (pnInter * ratioPs)

	if deleterious > pn
		deleterious = 0
		pnNeutral = round(pn - deleterious,digits=3)
	else
		pnNeutral = round(pn - deleterious,digits=3)
	end
	# output['alpha'] = 1 - (((pn - deleterious) / ps) * (ds / dn))
	output["alpha"] = round(1 - ((pnNeutral/ps) * (ds / dn)),digits=3)
	#  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps),Int(ceil(pnNeutral)),Int(ds),Int(dn)))


	if (!isnothing(m))

		mn = m[1]; ms = m[2]
		# ## Estimation of b: weakly deleterious
		output["b"] = (deleterious / ps) * (ms / mn)

		## Estimation of f: neutral sites
		output["f"] = (ms * pnNeutral) / (mn * ps)

		## Estimation of d, strongly deleterious sites
		output["d"] = 1 - (output["f"] + output["b"])

		ka      = dn / mn
		ks       = ds / ms
		output["omega"]    = ka/ ks


		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end

"""
	fwwMK(sfs,divergence,m,cutoff)

Function to estimate the asmpytotic value of α(x).

# Arguments
 - `alphaValues::Array{Float64,1}`: α(x) array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function fwwMK(;sfs::Array,divergence::Array,m::T,cutoff::Float64=0.15) where {T<:Union{Nothing,Array}}

	output = OrderedDict{String,Float64}()

	ps = sum(sfs[:,2])
	pn = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	
	deleterious = 0
	### Estimating slightly deleterious with pn/ps ratio
	fltInter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
	pnInter = sum(sfs[fltInter,2])
	psInter = sum(sfs[fltInter,3])


	# output['alpha'] = 1 - (((pn - deleterious) / ps) * (ds / dn))
	output["alpha"] = round(1 - ((pnInter/psInter) * (ds / dn)),digits=3)
	#  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(psInter),Int(ceil(pnInter)),Int(ds),Int(dn)))


	if (!isnothing(m))

		mn = m[1]; ms = m[2]
		ka      = dn / mn
		ks       = ds / ms
		output["omega"]    = ka/ ks


		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end

"""
	standardMK(x,columns)

Function to estimate the original α value

# Arguments
 - `alphaValues::Array{Float64,1}`: α(x) array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function standardMK(;sfs::Array,divergence::Array,m::T=nothing) where {T<:Union{Nothing,Array}}

	output = OrderedDict{String,Float64}()

	pn = sum(sfs[:,2])
	ps = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	
	output["alpha"] = round(1 - ((pn/ps) * (ds / dn)),digits=3)
	#  method = :mnnlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps),Int(ceil(pn)),Int(ds),Int(dn)))

	if (!isnothing(m))

		mn = m[1]; ms = m[2]

		ka      = dn / mn
		ks       = ds / ms
		output["omega"]    = ka/ ks


		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end
