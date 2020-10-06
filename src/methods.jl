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

	return [asymp ciLow ciHigh]
	# return asymp
end
