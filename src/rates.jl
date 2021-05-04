"""
	rates(param::parameters,iterations::Int64,divergence::Array,sfs::Array)

Function to solve randomly *N* scenarios. N = iterations ⋅ param.bRange

# Arguments
 - `param::parameters`
 - `convolutedSamples::binomialDict`
 - `gH::Array{Int64,1}`
 - `gL::Union{Array{Int64,1},Nothing}`
 - `gamNeg::Array{Int64,1}`
 - `shape::Float64=0.184`
 - `iterations::Int64`
 - `output::String`
# Returns
 - `Array`: summary statistics.
 - `Output`: HDF5 file containing models solved and rates.
"""
function rates(;param::parameters,convolutedSamples::binomialDict,gH::Array{Int64,1},gL::Union{Array{Int64,1},Nothing},gamNeg::Array{Int64,1},theta::Union{Float64,Nothing}=nothing,rho::Union{Float64,Nothing}=nothing,shape::Float64=0.184,iterations::Int64,output::String)

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

	# Creating N models to iter in threads. Set N models (paramerters) and sampling probabilites (binomialDict)
	nParam  = [param for i in 1:iterations];
	nBinom  = [convolutedSamples for i in 1:iterations];
	
	# Random strong selection coefficients
	ngh     = rand(repeat(gH,iterations),iterations);
	# Random negative selection coefficients
	ngamNeg = rand(repeat(gamNeg,iterations),iterations);
	
	# Random θ on coding regions
	if !isnothing(theta)
		θ = fill(theta,iterations)
	else
		θ = rand(0.0005:0.0005:0.01,iterations)
	end
	# Random ρ on coding regions
	if !isnothing(rho)
		ρ = fill(rho,iterations)
	else
		ρ = rand(0.0005:0.0005:0.05,iterations)
	end
	
	# Estimations to thread pool. 
	# Allocate ouput in SharedArray
	out    = SharedArray{Float64,3}(size(param.bRange,2),(size(param.dac,1) *2) + 12,iterations)
	@time @sync @distributed for i in 1:iterations
		# Each iteration solve 1 model accounting all B value in param.bRange
		@inbounds out[:,:,i] = iterRates(nParam[i], nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i]);
	end

	# Remove the workers to free memory resources
	# SharedArray is not remove after this process
	for i in Distributed.workers()
		rmprocs(i)
	end

	# Reducing array
	df = vcat(eachslice(out,dims=3)...);
	
	# Saving models and rates
	models = DataFrame(df[:,1:8],[:B,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:ρ])
	neut   = df[:,9:(8+size(param.dac,1))]
	sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)]
	dsdn   = df[:,(end-3):end]

	# Saving multiple summary statistics
	n = OrderedDict{Int,Array}()
	s = OrderedDict{Int,Array}()
	for i in eachindex(param.dac)
		n[param.dac[i]] = neut[:,i]
		s[param.dac[i]] = sel[:,i]
	end

	# Writting HDF5 file
	JLD2.jldopen(output, "a+") do file
		file[string(param.N)* "/" * string(param.n) * "/models"] = models
		file[string(param.N)* "/" * string(param.n) * "/neut"]   = n
		file[string(param.N)* "/" * string(param.n) * "/sel"]    = s
		file[string(param.N)* "/" * string(param.n) * "/dsdn"]   = dsdn
		file[string(param.N)* "/" * string(param.n) * "/dac"]    = param.dac
	end

	return df
end

"""
	iterRates(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,divergence::Array,sfs::Array)

Estimating rates given a model for all B range.

# Arguments
 - `param::parameters`
 - `convolutedSamples::binomialDict`
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
function iterRates(param::parameters,convolutedSamples::binomialDict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64)

	# Creating model to solve
	# Γ distribution
	param.al    = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
	# α, αW
	param.alLow = alLow; param.alTot = alTot;
	# Positive selection coefficients
	param.gH    = gH;param.gL = gL
	# Mutation rate and recomb
	param.thetaMidNeutral = θ; param.rho = ρ
	# Solving θ on non-coding region and probabilites to get α value without BGS
	param.B = 0.999
	setThetaF!(param)
	setPpos!(param)

	# Allocate array to solve the model for all B values
	r = spzeros(size(param.bRange,2),(size(param.dac,1) * 2) + 12)
	for j in eachindex(param.bRange)
		# Set B value
		param.B = param.bRange[j]
		# Solve θ non-coding for the B value.
		setThetaF!(param)
		# Solve model for the B value
		@inbounds r[j,:] = gettingRates(param,convolutedSamples.bn[param.B])
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
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN),false)
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	## Outputs
	#=α = @. 1 - (ds/dn) * (sel/neut)=#

	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	#=fN_nopos       = fN*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fNeg_nopos     = fNeg*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fPosL_nopos    = fPosL*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fPosH_nopos    = fPosH*(param.thetaMidNeutral/2.)*param.TE*param.NN

	ds_nopos       = fN_nopos
	dn_nopos       = fNeg_nopos + fPosL_nopos + fPosH_nopos
	dnS_nopos      = dn_nopos - fPosL_nopos

	## Polymorphism
	sel_nopos = selN

	## Outputs
	αW         = param.alLow/param.alTot
	α_nopos    = @. 1 - (ds_nopos/dn_nopos) * (sel_nopos/neut)=#


	#=alphas = round.(vcat(α_nopos[param.dac[end]] * αW , α_nopos[param.dac[end]] * (1 - αW), α_nopos[param.dac[end]]), digits=5)=#

	##########
	# Output #
	##########
	#=analyticalValues::Array{Float64,2} = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH)'=#
	analyticalValues::Array{Float64,2} = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH)'

	return (analyticalValues)
end