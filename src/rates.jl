"""
	rates(param::parameters,iterations::Int64,divergence::Array,sfs::Array)

Function to solve randomly *N* scenarios. The function will create *N* models, defined by ```Analytical.parameters()```, to estimate analytically fixation and polymorphic rates for each model. The rates will be used to compute summary statistics required at ABC. The function output a HDF5 file containing the solved models, the selected DAC and the analytical rates. 

If rho and/or theta are set to ```nothing```, the function will input random values given the range 0.0005:0.0005:0.01. Otherwise you can fix the values.

If gL is set to ```nothing```, the function will not account the role of the weakly selected alleles in the estimation.

# Arguments
 - `param::parameters`: mutable structure containing the model
 - `convoluted_samples::binomial_dict` : structure containing the binomial convolution
 - `gH::Array{Int64,1}` : Range of strong selection coefficients
 - `gL::Union{Array{Int64,1},Nothing}`: Range of weak selection coefficients
 - `gamNeg::Array{Int64,1}` : Range of deleterious selection coefficients
 - `theta::Union{Float64,Nothing}` : Population-scaled mutation rate on coding region
 - `rho::Union{Float64,Nothing}` : Population-scaled recombination rate
 - `shape::Float64=0.184` : DFE shape parameter
 - `iterations::Int64` : Number of solutions
 - `output::String` : File to output HDF5 file
# Returns
 - `Array`: summary statistics
 - `Output`: HDF5 file containing models solved and rates.
"""
function rates(;param::parameters,
				convoluted_samples::binomial_dict,
				gH::S,
				gL::S,
				gamNeg::S,
				theta::Union{Float64,Nothing}=0.001,
				rho::Union{Float64,Nothing}=0.001,
				shape::Float64=0.184,
				iterations::Int64,
				output::String,
				scheduler::String="local") where S <: Union{Array{Int64,1},UnitRange{Int64},Nothing}
	
	# Iterations = models to solve
	# Factor to modify input Γ(shape) parameter. Flexible Γ distribution over negative alleles
	fac     = rand(-2:0.05:2,iterations)
	afac    = @. param.al*(2^fac)
	
	# Deleting shape > 1. Negative alpha_x values
	idx = findall(afac .> 1)
	if !isempty(idx)
		afac[idx] = rand(afac[afac .< 1],size(idx,1))
	end

	# Random α values
	nTot    = rand(0.1:0.01:0.9,iterations)
	
	# Defining αW. It is possible to solve non-accounting for weak fixations
	if isnothing(gL)
		# Setting αW to 0 for all estimations
		nLow    = fill(0.0,iterations)
		# Random strong selection coefficients
		ngl     = rand(repeat([1],iterations),iterations);
	else
		# Setting αW as proportion of α
		lfac    = rand(0.0:0.05:0.9,iterations)
		nLow    = @. nTot * lfac
		# Random weak selection coefficients
		ngl     = rand(repeat(gL,iterations),iterations);
	end

	# Creating N models to iter in threads. Set N models (paramerters) and sampling probabilites (binomial_dict)
	nParam  = [param for i in 1:iterations];
	nBinom  = [convoluted_samples for i in 1:iterations];
	
	# Random strong selection coefficients
	ngh     = rand(repeat(gH,iterations),iterations);
	# Random negative selection coefficients
	ngamNeg = rand(repeat(gamNeg,iterations),iterations);

	# Random θ on coding regions
    #=obsNeut        = sfs.p0./sum(sfs.p0) 
    n              = 1 ./collect(1:(param.nn-1))
    expNeut        = n ./ sum(n)=#

	if !isnothing(theta)
		θ = fill(theta,iterations)
		#=θᵣ= fill(theta .* obsNeut ./ expNeut,iterations)=#
	else
		θ = rand(0.0005:0.0005:0.01,iterations)
		#=θᵣ= map(x -> x .* obsNeut ./ expNeut,θ)=#
	end

	# Random ρ on coding regions
	if !isnothing(rho)
		ρ = fill(rho,iterations)
	else
		ρ = rand(0.0005:0.0005:0.05,iterations)
	end
	
	# Estimations to distributed workers
	if scheduler == "local"
		println(scheduler)
		@time out = Distributed.pmap(iterRates,nParam, nBinom, nTot, nLow, ngh, ngl, ngamNeg, afac, θ, ρ);
		#=@time out = Distributed.pmap(iterRates,nParam, nBinom, nTot, nLow, ngh, ngl, ngamNeg, afac, θ, θᵣ, ρ);=#
	else
		println(scheduler)
		@time out = ParallelUtilities.pmapbatch( iterRates,nParam, nBinom, nTot, nLow, ngh, ngl, ngamNeg, afac, θ, ρ);
		#=@time out = ParallelUtilities.pmapbatch( iterRates,nParam, nBinom, nTot, nLow, ngh, ngl, ngamNeg, afac, θ, θᵣ, ρ);=#
	end

	# Remove the workers to free memory resources
	#=for i in workers()
		rmprocs(i)
	end=#

	# Reducing output array
	df = vcat(out...)
	
	# Saving models and rates
	models = DataFrame(df[:,1:8],[:B,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:ρ])
	neut   = df[:,9:(8+size(param.dac,1))]
	sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)]
	dsdn   = Array(df[:,(end-5):end-2])

	sum_pol = df[:,end-1:end]

	# Saving multiple summary statistics
	n = OrderedDict{Int,Array}()
	s = OrderedDict{Int,Array}()
	for i in eachindex(param.dac)
		n[param.dac[i]] = neut[:,i]
		s[param.dac[i]] = sel[:,i]
	end

	# Writting HDF5 file
	JLD2.jldopen(output, "a+") do file
		file[string(param.N)* "/" * string(param.n) * "/models"] = models;
		file[string(param.N)* "/" * string(param.n) * "/neut"]   = n;
		file[string(param.N)* "/" * string(param.n) * "/sel"]    = s;
		file[string(param.N)* "/" * string(param.n) * "/dsdn"]   = dsdn;
		file[string(param.N)* "/" * string(param.n) * "/pol"]   = sum_pol;
		file[string(param.N)* "/" * string(param.n) * "/dac"]    = param.dac;
	end;
end

"""
	iterRates(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,divergence::Array,sfs::Array)

Estimating rates given a model for all B range.

# Arguments
 - `param::parameters`
 - `convoluted_samples::binomial_dict`
 - `alTot::Float64`
 - `alLow::Float64`
 - `gH::Int64`
 - `gL::Int64`
 - `gamNeg::Int64`
 - `afac::Float64`
 - `ρ::Float64`
 - `θ::Float64`
# Output
 - `Array{Float64,2}`
"""
function iterRates(param::parameters,convoluted_samples::binomial_dict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64)

	# Creating model to solve
	# Γ distribution
	param.al    = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
	# α, αW
	param.alLow = alLow; param.alTot = alTot;
	# Positive selection coefficients
	param.gH    = gH;param.gL = gL
	# Mutation rate and recomb
	param.thetaMidNeutral = θ; param.rho = ρ
	#=param.thetaMidNeutral = θ; param.θᵣ .= θᵣ; param.rho = ρ=#
	# Solving θ on non-coding region and probabilites to get α value without BGS
	param.B = 0.999
	setThetaF!(param)
	setPpos!(param)

	# Allocate array to solve the model for all B values
	r = zeros(size(param.B_bins,1),(size(param.dac,1) * 2) + 14)
	for j in eachindex(param.B_bins)
		# Set B value
		param.B = param.B_bins[j]
		# Solve θ non-coding for the B value.
		setThetaF!(param)
		# Solve model for the B value
		tmp = try
			gettingRates(param,convoluted_samples.bn[param.B])
		catch e
			zeros(size(param.dac,1) * 2 + 14)'
		end
		@inbounds r[j,:] = tmp
	end

	return r
end

"""
	gettingRates(gammaL,gammaH,pposL,pposH,observedData,nopos)

Estimating analytical rates of fixation and polymorphism to approach α value accouting for background selection, weakly and strong positive selection. Output values will be used to sample from a Poisson distribution the total counts of polymorphism and divergence using observed data. 

# Arguments
 - `param::parameters`
 - `cnvBinom::SparseMatrixCSC{Float64,Int64}`
# Returns
 - `Array{Float64,2}` containing solved model, fixation and polymorphic rates
"""
function gettingRates(param::parameters,cnvBinom::SparseMatrixCSC{Float64,Int64})

	################################################
	# Subset rates accounting for positive alleles #
	################################################

	# Fixation
	fN       = param.B*fixNeut(param)
	fNeg     = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
	fPosL    = fixPosSim(param,param.gL,0.5*param.pposL)
	fPosH    = fixPosSim(param,param.gH,0.5*param.pposH)

	ds       = fN
	dn       = fNeg + fPosL + fPosH

	# Polymorphism
	neut::Array{Float64,1} = DiscSFSNeutDown(param,cnvBinom)
	selH::Array{Float64,1} = if isinf(exp(param.gH * 2))
		DiscSFSSelPosDownArb(param,param.gH,param.pposH,cnvBinom)
	else
		DiscSFSSelPosDown(param,param.gH,param.pposH,cnvBinom)
	end
	selL::Array{Float64,1} = DiscSFSSelPosDown(param,param.gL,param.pposL,cnvBinom)
	selN::Array{Float64,1} = DiscSFSSelNegDown(param,param.pposH+param.pposL,cnvBinom)
	
	# Cumulative rates
	tmp = cumulative_sfs(hcat(neut,selH,selL,selN),false)
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	##########
	# Output #
	##########
	analyticalValues::Array{Float64,2} = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH,neut[1],sel[1])'

	return (analyticalValues)
end


###########################
function rates_threads(;param::parameters,
                convoluted_samples::binomial_dict,
                gH::S,
                gL::S,
                gamNeg::S,
                theta::Union{Float64,Nothing}=0.001,
                rho::Union{Float64,Nothing}=0.001,
                shape::Float64=0.184,
                iterations::Int64,
                output::String) where S <: Union{Array{Int64,1},UnitRange{Int64},Nothing}
    

    # Iterations = models to solve
    # Factor to modify input Γ(shape) parameter. Flexible Γ distribution over negative alleles
    fac     = rand(-2:0.05:2,iterations)
    afac    = @. param.al*(2^fac)
    
    # Deleting shape > 1. Negative alpha_x values
    idx = findall(afac .> 1)
    if !isempty(idx)
        afac[idx] = rand(afac[afac .< 1],size(idx,1))
    end

    # Random α values
    nTot    = rand(0.1:0.01:0.9,iterations)
    
    # Defining αW. It is possible to solve non-accounting for weak fixations
    if isnothing(gL)
        # Setting αW to 0 for all estimations
        nLow    = fill(0.0,iterations)
        # Random strong selection coefficients
        ngl     = rand(repeat([1],iterations),iterations);
    else
        # Setting αW as proportion of α
        lfac    = rand(0.0:0.05:0.9,iterations)
        nLow    = @. nTot * lfac
        # Random weak selection coefficients
        ngl     = rand(repeat(gL,iterations),iterations);
    end

    # Creating N models to iter in threads. Set N models (paramerters) and sampling probabilites (binomial_dict)
    nParam  = [param for i in 1:iterations];
    nBinom  = [convoluted_samples for i in 1:iterations];
    
    # Random strong selection coefficients
    ngh     = rand(repeat(gH,iterations),iterations);
    # Random negative selection coefficients
    ngamNeg = rand(repeat(gamNeg,iterations),iterations);

    # Random θ on coding regions
    #=obsNeut        = sfs.p0./sum(sfs.p0) 
    n              = 1 ./collect(1:(param.nn-1))
    expNeut        = n ./ sum(n)=#

    if !isnothing(theta)
        θ = fill(theta,iterations)
        #=θᵣ= fill(theta .* obsNeut ./ expNeut,iterations)=#
    else
        θ = rand(0.0005:0.0005:0.01,iterations)
        #=θᵣ= map(x -> x .* obsNeut ./ expNeut,θ)=#
    end

    # Random ρ on coding regions
    if !isnothing(rho)
        ρ = fill(rho,iterations)
    else
        ρ = rand(0.0005:0.0005:0.05,iterations)
    end

    ######Parallel; 10^5 -> 62.76"
    x = zeros(1,3,iterations);
    @inbounds @sync for i in 1:iterations
        Base.Threads.@spawn begin
            @inbounds x[:,:,i] = solve(param,convoluted_samples, nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i])
        end
    end    
    
    #y = ThreadsX.map(solve, nParam,nBinom, nTot, nLow, ngh, ngl, ngamNeg, afac, θ, ρ)

    x = repeat(x, outer = [size(param.B_bins,1), 1, 1]);
    x = vcat([ @view x[:,:,i] for i=1:iterations]...);
    b = repeat(param.B_bins,iterations);
    e = hcat(x,b);

    out = zeros(iterations,32,size(param.B_bins,1));

    @showprogress for (j,val) in enumerate(reverse(param.B_bins))
        tmp = e[e[:,4] .== val,:]
        for i in 1:iterations
            Base.Threads.@spawn begin
               @inbounds out[i,:,j] = r(nParam[i],nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i],tmp[i,:]);
            end
        end
    end

    df = vcat([@view out[:,:,i] for i=1:size(out,3)]...);

    # Saving models and rates
    models = DataFrame(df[:,1:8],[:B,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:ρ]);
    neut   = df[:,9:(8+size(param.dac,1))];
    sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)];
    dsdn   = Array(df[:,(end-5):end-2]);

    sum_pol = df[:,end-1:end];

    # Saving multiple summary statistics
    n = OrderedDict{Int,Array}();
    s = OrderedDict{Int,Array}();
    for i in eachindex(param.dac)
        n[param.dac[i]] = neut[:,i]
        s[param.dac[i]] = sel[:,i]
    end

    #=return models=#
    # Writting HDF5 file
    JLD2.jldopen(output, "a+") do file
        file[string(param.N)* "/" * string(param.n) * "/models"] = models;
        file[string(param.N)* "/" * string(param.n) * "/neut"]   = n;
        file[string(param.N)* "/" * string(param.n) * "/sel"]    = s;
        file[string(param.N)* "/" * string(param.n) * "/dsdn"]   = dsdn;
        file[string(param.N)* "/" * string(param.n) * "/pol"]   = sum_pol;
        file[string(param.N)* "/" * string(param.n) * "/dac"]    = param.dac;
    end;
end

function solve(param::parameters,convoluted_samples::binomial_dict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64)
        # Creating model to solve
    # Γ distribution
    param.al = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
    # α, αW
    param.alLow = alLow; param.alTot = alTot;
    # Positive selection coefficients
    param.gH    = gH;param.gL = gL
    # Mutation rate and recomb
    param.thetaMidNeutral = θ; param.rho = ρ
    #=0.001.thetaMidNeutral = θ; param.θᵣ .= θᵣ; param.rho = ρ=#
    # Solving θ on non-coding region and probabilites to get α value without BGS
    param.B = 0.999
    setThetaF!(param)
    setPpos!(param)    

    return([param.pposL,param.pposH,param.thetaF])

end

function r(param::parameters,convoluted_samples::binomial_dict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64,e::Vector{Float64})

    # Creating model to solve
    # Γ distribution
    param.al    = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
    # α, αW
    param.alLow = alLow; param.alTot = alTot;
    # Positive selection coefficients
    param.gH    = gH;param.gL = gL
    # Mutation rate and recomb
    param.thetaF = e[3]; param.thetaMidNeutral = θ; param.rho = ρ
    # Fixation probabilites
    param.pposL = e[1]; param.pposH = e[2]

    # Solving θ on non-coding region and probabilites to get α value without BGS
    param.B = e[4]
    cnvBinom = convoluted_samples.bn[e[4]]
    # Solve θ non-coding for the B value.
    setThetaF!(param)
    # Solve model for the B value

    # Fixation
    fN       = param.B*fixNeut(param)
    fNeg     = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
    fPosL    = fixPosSim(param,param.gL,0.5*param.pposL)
    fPosH    = fixPosSim(param,param.gH,0.5*param.pposH)

    ds       = fN
    dn       = fNeg + fPosL + fPosH

    # Polymorphism
    neut = DiscSFSNeutDown(param,cnvBinom)
    selH = if isinf(exp(param.gH * 2))
        DiscSFSSelPosDownArb(param,param.gH,param.pposH,cnvBinom)
    else
        DiscSFSSelPosDown(param,param.gH,param.pposH,cnvBinom)
    end
    selL = DiscSFSSelPosDown(param,param.gL,param.pposL,cnvBinom)
    selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,cnvBinom)

    # Cumulative rates
    tmp = cumulative_sfs(hcat(neut,selH,selL,selN),false)
    splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
    neut, selH, selL, selN = splitColumns(tmp)
    sel = (selH+selL)+selN

    ##########
    # Output #
    ##########
    r = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH,neut[1],sel[1])'
    
    return r
end