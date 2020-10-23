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
function summaryStats(;param::parameters,alpha::Float64,shape::Float64=0.184,scale::Float64=0.000402,divergence::Array,sfs::Array,bins::Int64,dac::Array{Int64,1},iterations::Int64)

	iterations  = trunc(Int,iterations/17) + 1
	# N random prior combinations
    fac         = rand(-2:0.05:2,iterations,2)
    afac        = @. shape*(2^fac[:,1])
    bfac        = @. scale*(2^fac[:,2])
    alTot       = rand(collect(0.1:0.05:alpha),iterations)
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
	df  = reduce(vcat,tmp)
    idx = vcat(1:3,3 .+ dac)
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
function bgsIter(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,divergence::Array,sfs::Array,bins::Int64)

	# Matrix and values to solve
	dm 			= size(divergence,1)
    r           = Array{Float64}(undef, 17 * dm , bins + 3)
    param.al    = afac; param.be = bfac;
    param.alLow = alLow; param.alTot = alTot;

    # Solve probabilites without B effect to achieve α value
    param.B = 0.999
    set_theta_f!(param)
    setPpos!(param)

	iter = 1
    for j in eachindex(param.bRange)
        param.B = param.bRange[j]
        # Solve mutation given a new B value.
        set_theta_f!(param)
        # Solven given same probabilites probabilites ≠ bgs mutation rate.
        x,y,z::Array{Float64,2} = alphaByFrequencies(param,divergence,sfs,bins,0.999)
		# push!(r,z)
        r[iter:(iter + (dm - 1)),:] = z
		iter = iter + dm
    end

    return r
end
