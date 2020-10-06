"""
	summStats(param::parameters,iterations::Int64,divergence::Array,sfs::Array)

Function to solve randomly *N* scenarios

# Arguments
 - `param::parameters`
 - `iterations::Int64`
 - `divergence::Array`
 - `sfs::Array`
# Returns
 - `Array`: summary statistics
"""
function summaryStats(;param::parameters,alpha::Float64,shape::Float64=0.184,scale::Float64=0.000402,divergence::Array{Int64,T},sfs::T,bins::Int64,iterations::Int64) where T <:Int64

	iterations  = trunc(Int,iterations/17) + 1
	# N random prior combinations
    fac         = rand(-2:0.05:2,iterations)
    afac        = @. shape*(2^fac)
    bfac        = @. scale*(2^fac)
    alTot       = rand(collect(0.01:0.01:alpha),iterations)
    lfac        = rand(collect(0.1:0.1:0.9),iterations)
    alLow       = @. round(alTot * lfac,digits=5)
    nParam      = [param for i in 1:iterations]
    ndivergence = [divergence for i in 1:iterations]
    nSfs        = [sfs for i in 1:iterations]
    nBins       = [bins for i in 1:iterations]

	# Estimations to thread pool
    wp  = Distributed.CachingPool(Distributed.workers())
    tmp = Distributed.pmap(bgsIter,wp,nParam,afac,bfac,alTot,alLow,ndivergence,nSfs,nBins);

	# Output
	df = reduce(vcat,tmp)

	return df
end

"""
	bgsIter(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,divergence::Array,sfs::Array)

Function to input and solve one scenario given *N* background selection values (*B*). It is used at *summStats* function to solve all the scenarios in multi-threading. Please use *Distributed* module to add *N* threads using *addprocs(N)*. In addition you must to declare our module in all the threads using the macro *@everywhere* before to start. Check the example if don't know how to do it.

# Arguments
 - `param::parameters`
 - `iterations::Int64`
 - `divergence::Array`
 - `sfs::Array`
# Returns
 - `Array`: summary statistics
"""
function bgsIter(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,divergence::Array{Float64,1},sfs::Array{Float64,1},bins::Int64)

	# Matrix and values to solve
    r           = Array{Float64}(undef, 17, bins + 3)
    param.al    = afac; param.be = bfac;
    param.alLow = alLow; param.alTot = alTot;

    # Solve probabilites without B effect to achieve α value
    param.B = 0.999
    Analytical.set_theta_f!(param)
    Analytical.setPpos!(param)

    iter = 1
    for j in param.bRange

        param.B = j
        # Solve mutation given a new B value.
        Analytical.set_theta_f!(param)
        # Solven given same probabilites probabilites ≠ bgs mutation rate.
        x,y,z = Analytical.alphaByFrequencies(param,divergence,sfs,100,0.9)

        r[iter,:] = z
        iter = iter + 1;
    end

    return r
end
