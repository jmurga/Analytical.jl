"""
	asympFit(x,columns)

Function to estimate the asmpytotic value of α(x).

# Arguments
 - `alphaValues::Array{Float64,1}`: α(x) array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function amk(α_x::Array{Float64,1})

	# Model
	asymp_fit(x,p) = @. p[1] + p[2]*exp(-x*p[3])

	# Fit values
	fitted1    = curve_fit(asymp_fit,collect(1:size(α_x,1)),α_x,[-1.0,-1.0,1.0];lower=[-1.0,-1.0,1.0],upper=[1.0, 1.0, 10.0])
	fitted2    = curve_fit(asymp_fit,collect(1:size(α_x,1)),α_x,fitted1.param)
	asymp      = asymp_fit(size(α_x,1),fitted2.param)

	ci_low, ci_high   = try
		confidence_interval(fitted2)[1][1],LsqFit.confidence_interval(fitted2)[1][2]
	catch err
		(0.0,0.0)
	end

	return asymp,[ci_low,ci_high],fitted2.param
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

	pn = sum(sfs[:,2])
	ps = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	deleterious = 0
	### Estimating slightly deleterious with pn/ps ratio
	flt_low        = (sfs[:, 1] .<= cutoff)
	pn_low         = sum(sfs[flt_low,2])
	ps_low         = sum(sfs[flt_low,3])

	flt_inter      = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
	pn_inter       = sum(sfs[flt_inter,2])
	ps_inter       = sum(sfs[flt_inter,3])

	ratio_ps       = ps_low / ps_inter
	deleterious    = pn_low - (pn_inter * ratio_ps)

	if deleterious > pn
		deleterious = 0
		pn_neutral  = round(pn - deleterious,digits=3)
	else
		pn_neutral = round(pn - deleterious,digits=3)
	end
	# output['alpha'] = 1 - (((pn - deleterious) / ps) * (ds / dn))
	output["alpha"]  = round(1 - ((pn_neutral/ps) * (ds / dn)),digits=3)
	#  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps),Int(ceil(pn_neutral)),Int(ds),Int(dn)))


	if (!isnothing(m))

		mn = m[1]; ms = m[2]
		# ## Estimation of b: weakly deleterious
		output["b"] = (deleterious / ps) * (ms / mn)
		## Estimation of f: neutral sites
		output["f"] = (ms * pn_neutral) / (mn * ps)

		## Estimation of d, strongly deleterious sites
		output["d"] = 1 - (output["f"] + output["b"])

		ka                 = dn / mn
		ks                 = ds / ms
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
	flt_inter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
	pn_inter  = sum(sfs[flt_inter,2])
	ps_inter  = sum(sfs[flt_inter,3])


	# output['alpha'] = 1 - (((pn - deleterious) / ps) * (ds / dn))
	output["alpha"]  = round(1 - ((pn_inter/ps_inter) * (ds / dn)),digits=3)
	#  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps_inter),Int(ceil(pn_inter)),Int(ds),Int(dn)))


	if (!isnothing(m))

		mn = m[1]; ms = m[2]
		ka                 = dn / mn
		ks                 = ds / ms
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
		ka              = dn / mn
		ks              = ds / ms
		output["omega"] = ka/ ks

		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end
