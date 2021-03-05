"""
	rates(param::parameters,iterations::Int64,divergence::Array,sfs::Array)

Function to solve randomly *N* scenarios

# Arguments
 - `param::parameters`
 - `iterations::Int64`
 - `divergence::Array`
 - `sfs::Array`
# Returns
 - `Array`: summary statistics
"""
function rates(;param::parameters,convolutedSamples::binomialDict,gH::Array{Int64,1},gL::Union{Array{Int64,1},Nothing},gamNeg::Array{Int64,1},shape::Float64=0.184,iterations::Int64,output::String)

	fac     = rand(-2:0.05:2,iterations)
	afac    = @. param.al*(2^fac)
	
	idx = findall(afac .> 1)
	
	if !isempty(idx)
		afac[idx] .= 0.999
	end

	nTot    = rand(0.1:0.01:0.9,iterations)
	
	if isnothing(gL)
		nLow    = fill(0.0,iterations)
		ngl     = rand(repeat([1],iterations),iterations);
	else
		lfac    = rand(0.0:0.05:0.9,iterations)
		nLow    = @. nTot * lfac
		ngl     = rand(repeat(gL,iterations),iterations);
	end

	nParam  = [param for i in 1:iterations];
	nBinom  = [convolutedSamples for i in 1:iterations];
	ngh     = rand(repeat(gH,iterations),iterations);
	ngamNeg = rand(repeat(gamNeg,iterations),iterations);
	
	# Estimations to thread pool
	out    = SharedArray{Float64,3}(size(param.bRange,2),(size(param.dac,1) *2) + 12,iterations)
	@sync @distributed for i in 1:iterations
		@inbounds out[:,:,i] = iterRates(nParam[i], nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i]);
	end
	df = vcat(eachslice(out,dims=3)...);
	
	models = DataFrame(df[:,1:8],[:B,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:be])
	neut   = df[:,9:(8+size(param.dac,1))]
	sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)]
	dsdn   = df[:,(end-3):end]

	JLD2.jldopen(output, "a+") do file
		file[string(param.N)* "/" * string(param.n) * "/models"] = models
		file[string(param.N)* "/" * string(param.n) * "/neut"]   = neut
		file[string(param.N)* "/" * string(param.n) * "/sel"]    = sel
		file[string(param.N)* "/" * string(param.n) * "/dsdn"]   = dsdn
		#=file[string(param.N)* "/" * string(param.n) * "/alphas"] = alphas=#
		file[string(param.N)* "/" * string(param.n) * "/dac"]    = param.dac
	end

	return df
end

"""
	iterRates(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,divergence::Array,sfs::Array)
"""
function iterRates(param::parameters,convolutedSamples::binomialDict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64)

	# Matrix and values to solve
	param.al    = afac; param.be = abs(afac/gamNeg);
	param.alLow = alLow; param.alTot = alTot;
	param.gH    = gH;param.gL = gL; param.gamNeg = gamNeg

	# Solve probabilites without B effect to achieve α value
	param.B = 0.999
	setThetaF!(param)
	setPpos!(param)

	r = zeros(size(param.bRange,2),(size(param.dac,1) * 2) + 12)
	for j in eachindex(param.bRange)
		param.B = param.bRange[j]
		# Solve mutation given a new B value.
		setThetaF!(param)
		# Solven given same probabilites probabilites ≠ bgs mutation rate.
		#x,y,z::Array{Float64,2} = alphaByFrequencies(param,divergence,sfs,dac)
		@inbounds r[j,:] = gettingRates(param,convolutedSamples.bn[param.B])
	end
	return r
end

"""
	gettingRates(gammaL,gammaH,pposL,pposH,observedData,nopos)

Analytical α(x) estimation. We used the expected rates of divergence and polymorphism to approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. We solve α(x) from empirical observed values. The values will be use to sample from a Poisson distribution the total counts of polymorphism and divergence using the rates. The mutation rate, the locus length and the time of the branch should be proportional to the observed values.

```math
\\mathbb{E}[\\alpha_{x}] =  1 - \\left(\\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{N}]}\\frac{\\mathbb{E}[P_{N}]}{\\mathbb{E}[P_{S}]}\\right)
```

# Arguments
 - `gammaL::Int64`: strength of weakly positive selection
 - `gammaH::Int64`: strength of strong positive selection
 - `pposL`::Float64: probability of weakly selected allele
 - `pposH`::Float64: probability of strong selected allele
 - `observedData::Array{Any,1}`: Array containing the total observed divergence, polymorphism and site frequency spectrum.
 - `nopos::String("pos","nopos","both")`: string to perform α(x) account or not for both positive selective alleles.

# Returns
 - `Tuple{Array{Float64,1},Array{Float64,2}}` containing α(x) and the summary statistics array (Ds,Dn,Ps,Pn,α).
"""
function gettingRates(param::parameters,cnvBinom::SparseMatrixCSC{Float64,Int64})

	##############################################################
	# Accounting for positive alleles segregating due to linkage #
	##############################################################

	# Fixation
	fN       = param.B*fixNeut(param)
	fNeg     = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
	fPosL    = fixPosSim(param,param.gL,0.5*param.pposL)
	fPosH    = fixPosSim(param,param.gH,0.5*param.pposH)

	ds       = fN
	dn       = fNeg + fPosL + fPosH

	## Polymorphism	## Polymorphism
	neut::Array{Float64,1} = DiscSFSNeutDown(param,cnvBinom)

	selH::Array{Float64,1} = if isinf(exp(param.gH * 2))
		DiscSFSSelPosDown(param,param.gH,param.pposH,cnvBinom)
	else
		DiscSFSSelPosDownArb(param,param.gH,param.pposH,cnvBinom)
	end

	selL::Array{Float64,1} = DiscSFSSelPosDown(param,param.gL,param.pposL,cnvBinom)
	selN::Array{Float64,1} = DiscSFSSelNegDown(param,param.pposH+param.pposL,cnvBinom)
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN),false)
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	## Outputs
	α = @. 1 - (ds/dn) * (sel/neut)

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

	##########
	# Output #
	##########

	#=alphas = round.(vcat(α_nopos[param.dac[end]] * αW , α_nopos[param.dac[end]] * (1 - αW), α_nopos[param.dac[end]]), digits=5)=#
	analyticalValues::Array{Float64,2} = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.be,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH)'

	return (analyticalValues)
end