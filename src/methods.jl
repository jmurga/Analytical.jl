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
	imputedMK(x,columns)

Function to estimate the asmpytotic value of α(x).

# Arguments
 - `alphaValues::Array{Float64,1}`: α(x) array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function imputedMK(;sfs::Array{Float64,2},divergence::T,m::Union{Nothing,T}=nothing,cutoff::Float64=0.15) where {T<:Union{Array{Float64,1},Array{Int64,1}}}

    output = OrderedDict{String,Float64}()

    pi = sum(sfs[:,2])
    p0 = sum(sfs[:,3])
    di = divergence[1]
    d0 = divergence[2]
    
    
    deleterious = 0
    ### Estimating slightly deleterious with pi/p0 ratio
    fltLow = (sfs[:, 1] .<= cutoff)
    piLow   = sum(sfs[fltLow,2])
    p0Low   = sum(sfs[fltLow,3])

    fltInter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
    piInter = sum(sfs[fltInter,2])
    p0Inter = sum(sfs[fltInter,3])

    ratioP0       = p0Low / p0Inter
    deleterious   = piLow - (piInter * ratioP0)
    piNeutral     = round(pi - deleterious,digits=3)

    # output['alpha'] = 1 - (((pi - deleterious) / p0) * (d0 / di))
    output["alpha"] = round(1 - ((piNeutral/p0) * (d0 / di)),digits=3)
    #  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(p0),Int(ceil(piNeutral)),Int(d0),Int(di)))


    if (!isnothing(m))

    	mi = m[1]; m0 = m[2]
	    # ## Estimation of b: weakly deleterious
	    output["b"] = (deleterious / p0) * (m0 / mi)

	    ## Estimation of f: neutral sites
	    output["f"] = (m0 * piNeutral) / (mi * p0)

	    ## Estimation of d, strongly deleterious sites
	    output["d"] = 1 - (output["f"] + output["b"])

        # output["Ka"]       = di / mi
		# output["Ks"]       = d0 / m0
		# output["omega"]    = output["Ka"] / output["Ks"]


		## Omega A and Omega D
		# output["omegaA"] = output["omega"] * output["alpha"]
		# output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end
