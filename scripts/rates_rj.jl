import Parameters: @with_kw

@with_kw mutable struct parameters
	gamNeg::Int64              = -457
	gL::Int64                  = 10
	gH::Int64                  = 500
	alLow::Float64             = 0.2
	alTot::Float64             = 0.4
	thetaF::Float64            = 1e-3
	al::Float64                = 0.184
	be::Float64                = 0.000402
	B::Float64                 = 0.999 
	B_bins::Array{Float64,2}   = permutedims(push!(collect(0.1:0.025:0.975),0.999))
	pposL::Float64             = 0
	pposH::Float64             = 0.001
	N::Int64                   = 1000
	n::Int64                   = 500
	Lf::Int64                  = 2*10^5
	rho::Float64               = 0.001
	TE::Float64                = 5.0
	diploid::Bool              = false

	NN::Int64 = 2*N
	nn::Int64 = 2*n
	thetaMidNeutral::Array{Float64,1}   = fill(0.001,nn-1)
	dac::Array{Int64,1} = [1,2,4,5,10,20,50,200,500,700]

end


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
    obsNeut        = sfs.p0./sum(sfs.p0) 
    n              = 1 ./collect(1:(param.nn-1))
    expNeut        = n ./ sum(n)

	if !isnothing(theta)
		θ = fill(param.thetaMidNeutral .* obsNeut ./ expNeut,iterations)
	else
		θrng = rand(0.0005:0.0005:0.01,iterations)
		θrng = fill.(θrng,param.nn-1)
		θ = map(x -> x .* obsNeut ./ expNeut,θrng)
	end
	
	# Random ρ on coding regions
	if !isnothing(rho)
		ρ = fill(rho,iterations)
	else
		ρ = rand(0.0005:0.0005:0.05,iterations)
	end
	
	# Estimations to thread pool. 
	# Allocate ouput Array
	#=out    = SharedArray{Float64,3}(size(param.B_bins,2),(size(param.dac,1) *2) + 12,iterations)=#
	out    = SharedArray{Float64,3}(size(param.B_bins,2),(param.nn * 2 - 2) + 12,iterations)
	@time for i in 1:iterations
		# Each iteration solve 1 model accounting all B value in param.B_bins
		@inbounds out[:,:,i] = iterRates(nParam[i], nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i]);
	end
	#=param = nParam[i]; convolutedSamples=nBinom[i];alTot= nTot[i];alLow= nLow[i];gH= ngh[i];gL= ngl[i];gamNeg= ngamNeg[i];afac= afac[i];θ= θ[i]; ρ= ρ[i];obsNeut = sNeut[i]=#
	# Reducing array
	df = vcat(eachslice(out,dims=3)...);
	
	# Saving models and rates
	models = DataFrame(df[:,1:8],[:B,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:θ])
	neut   = df[:,9:param.nn-1+8]
	sel    = df[:,(9+param.nn-1):(8+param.nn*2-2)]
	dsdn   = df[:,(end-3):end]

	n = OrderedDict{Int,Array}()
	s = OrderedDict{Int,Array}()
	for i in param.dac
		n[i] = neut[:,i]
		s[i] = sel[:,i]
	end

	#=neut   = df[:,9:(8+size(param.dac,1))]
	sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)]
	dsdn   = df[:,(end-3):end]=#

	# Writting HDF5 file
	JLD2.jldopen(output, "a+") do file
		file[string(param.N)* "/" * string(param.n) * "/models"] = models
		file[string(param.N)* "/" * string(param.n) * "/neut"]   = neut
		file[string(param.N)* "/" * string(param.n) * "/sel"]    = sel
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
function iterRates(param::parameters,convolutedSamples::binomialDict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64,obsNeut::Array{Float64})

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
	#=r = zeros(size(param.B_bins,2),(size(param.dac,1) * 2) + 12)=#
	r = spzeros(size(param.B_bins,2),(param.nn * 2 - 2) + 12)
	for j in eachindex(param.B_bins)
		# Set B value
		param.B = param.B_bins[j]
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
	neut = DiscSFSNeutDown(param,cnvBinom)
	selH = if isinf(exp(param.gH * 2))
		DiscSFSSelPosDownArb(param,param.gH,param.pposH,cnvBinom)
	else
		DiscSFSSelPosDown(param,param.gH,param.pposH,cnvBinom)
	end
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL,cnvBinom)
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,cnvBinom)
	# Cumulative rates
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN),false)
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	##########
	# Output #
	##########

	analyticalValues::Array{Float64,2} = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut,sel,ds,dn,fPosL,fPosH)'

	return (analyticalValues)
end


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

Expected rate of neutral allele frequency reduce by backgrou	nd selection. The spectrum depends on the number of individual []

```math
\\mathbb{E}[Ps_{(x)}] = \\sum{x^{*}=x}{x^{*}=1}f_{B}(x)
```

# Return:
 - `Array{Float64}`: expected rate of neutral alleles frequencies.
"""
function DiscSFSNeutDown(param::parameters,binom::SparseMatrixCSC{Float64,Int64})

	NN2 = convert(Int64,ceil(param.NN*param.B))
	# Allocating variables

	neutralSfs(i::Int64) = 1.0/(i)

	x = collect(0:NN2)
	solvedNeutralSfs = neutralSfs.(x)
	replace!(solvedNeutralSfs, Inf => 0.0)

	# subsetDict = get(param.bn,param.B,1)
	# subsetDict = binom
	out::Array{Float64,1} = param.B.*(param.thetaMidNeutral).*0.25.*(binom*solvedNeutralSfs)
	# out = @view out[2:end-1]
	# out = out[2:end-1]

	return out
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
function DiscSFSSelPosDown(param::parameters,gammaValue::Int64,ppos::Float64,binom::SparseMatrixCSC{Float64,Int64})

	if ppos == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		redPlus = phiReduction(param,gammaValue)

		# Solving sfs
		NN2 = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*param.B

		gammaExp1 = exp(gammaCorrected*2)
		gammaExp2 = exp(gammaCorrected*-2)

		positiveSfs(i::Float64,g1::Float64=gammaExp1,g2::Float64=gammaExp2,ppos::Float64=ppos) = Float64(ppos*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i))))

		# Original
		# ppos*0.5*(ℯ^(2*gammaCorrected)*(1-ℯ^(-2.0*gammaCorrected*(1.0-i)))/((ℯ^(2*gammaCorrected)-1.0)*i*(1.0-i)))

		# Allocating outputs
		solvedPositiveSfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2))
		replace!(solvedPositiveSfs, NaN => 0.0)

		# subsetDict = get(param.bn,param.B,1)
		# out               = (param.thetaMidNeutral)*redPlus*0.75*(subsetDict*solvedPositiveSfs)
		out::Array{Float64,1} = (param.thetaMidNeutral).*redPlus.*0.75.*(binom*solvedPositiveSfs)
		# out = out[2:end-1]

	end

	return out
end

function DiscSFSSelPosDownArb(param::parameters,gammaValue::Int64,ppos::Float64,binom::SparseMatrixCSC{Float64,Int64})

	if ppos == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else
		redPlus = phiReduction(param,gammaValue)

		# Solving sfs
		NN2 = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*param.B
		gammaExp1::Arb = exp(Arb(gammaCorrected*2,prec=10))
		gammaExp2::Arb = exp(Arb(gammaCorrected*-2,prec=10))

		positiveSfs(i::Float64,g1::Arb=gammaExp1,g2::Arb=gammaExp2,ppos::Float64=ppos) = Float64(ppos*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i))))
		# Allocating outputs
		solvedPositiveSfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2))
		replace!(solvedPositiveSfs, NaN => 0.0)
		out::Array{Float64,1} = (param.thetaMidNeutral).*redPlus.*0.75.*(binom*solvedPositiveSfs)
		# out = out[2:end-1]

	end

	return out
end

######Slightly deleterious######
"""

	DiscSFSSelNegDown(param,ppos)

Expected rate of positive selected allele frequency reduce by background selection. Spectrum drawn on a gamma DFE. It depends on the number of individuals.

# Arguments
 - `ppos::Float64`: positive selected alleles probabilty.
# Return:
 - `Array{Float64}`: expected negative selected alleles frequencies.
"""
function DiscSFSSelNegDown(param::parameters,ppos::Float64,binom::SparseMatrixCSC{Float64,Int64})
	# subsetDict = get(param.bn,param.B,1)
	solvedNegative = DiscSFSSelNeg(param,ppos)
	# out::Array = param.B*(param.thetaMidNeutral)*0.75*(subsetDict*solvedNegative)
	out = param.B.*(param.thetaMidNeutral).*0.75.*(binom*solvedNegative)
	# out = @view out[2:end-1]

	# return out[2:end-1]
	return out
end

function DiscSFSSelNeg(param::parameters,ppos::Float64)

	beta     = param.be/(1.0*param.B)
	NN2      = convert(Int64, ceil(param.NN*param.B))
	xa       = collect(0:NN2)/NN2

	solveZ   = similar(xa)

	z(x::Float64,p::Float64=ppos) = (1.0-p)*(2.0^-param.al)*(beta^param.al)*(-SpecialFunctions.zeta(param.al,x+beta/2.0) + SpecialFunctions.zeta(param.al,(2+beta)/2.0))/((-1.0+x)*x)

	solveZ   = z.(xa)

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
function cumulativeSfs(sfsTemp::Array,freqs::Bool=true)

	out      = Array{Float64}(undef, size(sfsTemp,1),size(sfsTemp,2))

	if freqs
		idx = 2
	else
		idx = 1
	end

	out[1,idx:end] = sum(sfsTemp[:,idx:end],dims=1)

	@simd for i in 2:(size(sfsTemp)[1])

		#=app = view(out,i-1,:) .- view(sfsTemp,i-1,:)=#
		app = out[i-1,idx:end] .- sfsTemp[i-1,idx:end]

		if sum(app) > 0.0
			out[i,idx:end] = app
		else
			out[i,idx:end] = zeros(length(app))
		end
	end

	if freqs
		out[:,1] = sfsTemp[:,1]
	end
	
	return out
end

"""
	reduceSfs(sfsTemp,bins)

Function to reduce the SFS into N bins.
"""
function reduceSfs(sfsTemp::Array,bins::Int64)

	freq  = collect(0:(size(sfsTemp,1)-1))/size(sfsTemp,1)
	h1    = fit(Histogram,freq,0:(1/(bins-1)):1)
	xmap1 = StatsBase.binindex.(Ref(h1), freq)

	tmp = hcat(sfsTemp,xmap1)
	out = Array{Float64}(undef,bins-1 ,size(sfsTemp,2))
	vIter =  convert(Array,unique(xmap1)')
	@simd for i = eachindex(vIter)
		@inbounds out[i,2:end] = sum(tmp[tmp[:,end] .== i,2:end-1],dims=1)
	end

	out[:,1] = collect(1:(bins-1)) ./ bins

	return (out)
end

