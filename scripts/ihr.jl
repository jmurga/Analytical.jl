using Parameters, NLsolve, SpecialFunctions, Distributions, Roots, ArbNumerics, StatsBase, LsqFit, PoissonRandom, SparseArrays, Distributed, CSV, SharedArrays, JLD2, DataFrames, ProgressMeter

import GZip: open
import Parsers: parse
import OrderedCollections: OrderedDict
import FastaIO: readfasta

PATH = "/home/jmurga/.julia/dev/Analytical/src/"
include(PATH * "parameters.jl")
include(PATH * "fixations.jl")
include(PATH * "polymorphism.jl")
include(PATH * "summaryStatistics.jl")
include(PATH * "rates.jl")
include(PATH * "inferTools.jl")
include(PATH * "readFasta.jl")
include(PATH * "methods.jl")

param = parameters()
convolutedSamples = binomialDict()
binomOp!(param,convolutedSamples.bn)

function gettingRates(param::parameters,convolutedSamples::binomialDict)

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
	neut::Array{Float64,1} = DiscSFSNeutDown(param,convolutedSamples.bn[param.B])
	# neut = param.neut[param.B]

	selH::Array{Float64,1} = DiscSFSSelPosDown(param,param.gH,param.pposH,convolutedSamples.bn[param.B])
	selL::Array{Float64,1} = DiscSFSSelPosDown(param,param.gL,param.pposL,convolutedSamples.bn[param.B])
	selN::Array{Float64,1} = DiscSFSSelNegDown(param,param.pposH+param.pposL,convolutedSamples.bn[param.B])
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
	analyticalValues::Array{Float64,2} = cat(param.B,(fPosL+fPosH),param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.be,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH,dims=1)'

	return (analyticalValues)
end

function iterRates(;param::parameters,convolutedSamples::binomialDict,alTot::Float64,alLow::Float64,gH::Int64,gL=Int64,gamNeg::Int64,afac::Float64)

	# Matrix and values to solve
	dm 			= 1
	param.al    = afac; param.be = abs(afac/gamNeg);
	param.alLow = alLow; param.alTot = alTot;
	param.gH    = gH;param.gL = gL; param.gamNeg = gamNeg

	# Solve probabilites without B effect to achieve α value
	param.B = 0.999
	setThetaF!(param)
	setPpos!(param)

	r = zeros(size(param.bRange,2) * dm,(size(param.dac,1) * 2) + 13)
	for j in eachindex(param.bRange)
		param.B = param.bRange[j]
		# Solve mutation given a new B value.
		setThetaF!(param)
		# Solven given same probabilites probabilites ≠ bgs mutation rate.
		#x,y,z::Array{Float64,2} = alphaByFrequencies(param,divergence,sfs,dac)
		@inbounds r[j,:] = gettingRates(param,convolutedSamples)
	end
	r[:,2] = @. r[:,2]/r[end,2]
	return r
end


function rates(;param::parameters,convolutedSamples::binomialDict,gH::Array{Int64,1},gL::Array{Int64,1},gamNeg::Array{Int64,1},shape::Float64=0.184,iterations::Int64,output::String)

	#=iterations=1;gamNeg=[-500];gL=[5];gH=[500];shape=0.184=#
	fac    = rand(-2:0.05:2,iterations,2)
	afac   = @. param.al*(2^fac[:,1])
	#=bfac   = @. (param.al/param.be)*(2^fac[:,2])=#

	lfac   = rand(0.05:0.05:0.9,iterations)
	nTot   = rand(0.1:0.01:0.9,iterations)

	nLow   = @. nTot * lfac
	nParam = [param for i in 1:iterations];
	nBinom = [convolutedSamples for i in 1:iterations];
	ngh    = rand(repeat(gH,iterations),iterations);
	ngl    = rand(repeat(gL,iterations),iterations);
	ngamNeg    = rand(repeat(gamNeg,iterations),iterations);
	
	# Estimations to thread pool
	out    = SharedArray{Float64,3}(size(param.bRange,2),(size(param.dac,1) *2) + 13,iterations)
	@sync @distributed for i in eachindex(afac)
		@inbounds out[:,:,i] = iterRates(param = nParam[i],convolutedSamples=nBinom[i],alTot = nTot[i], alLow = nLow[i],gH=ngh[i],gL=ngl[i],gamNeg=ngamNeg[i],afac=afac[i]);
	end

	df = vcat(eachslice(out,dims=3)...);
	models = DataFrame(df[:,1:9],[:B,:hri,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:be])

    neut   = df[:,10:9+size(param.dac,1)]
    sel    = df[:,10+size(param.dac,1):9+size(param.dac,1)*2]
    dsdn   = df[:,end-3:end]

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


function summaryStatsFromRates(;param::parameters,rates::JLD2.JLDFile,divergence::Array,sfs::Array,summstatSize::Int64,replicas::Int64)

	tmp    = rates[string(param.N) * "/" * string(param.n)]
	idx    = StatsBase.sample.(fill(1:size(tmp["neut"],1),replicas),fill(summstatSize,replicas),replace=false)

	models = Array.(map(x -> view(tmp["models"],x,:), idx))
	neut   = Array.(map(x -> view(tmp["neut"],x,:), idx))
	sel    = Array.(map(x -> view(tmp["sel"],x,:), idx))
	neut   = Array.(map(x -> view(tmp["neut"],x,:), idx))
	dsdn   = Array.(map(x -> view(tmp["dsdn"],x,:), idx))

	expectedValues =  progress_pmap(samplingFromRates,models,sfs,divergence,neut,sel,dsdn);

	return(expectedValues)
end

function samplingFromRates(m::Array,s::Array,d::Array,nt::Array,sl::Array,x::Array)
	ds      = x[:,1]
	dn      = x[:,2]
	dweak   = x[:,3]
	dstrong = x[:,4]
	gn      = abs.(m[:,5])
	sh      = round.(m[:,end-1],digits=5)
	hri     = m[:,2]

	alxSummStat, alphasDiv, expectedDn, expectedDs, expectedPn, expectedPs = sampledAlpha(d=d,afs=s,λdiv=[ds,dn,dweak,dstrong],λpol=[permutedims(nt),permutedims(sl)])
	expectedValues = hcat(round.(alphasDiv,digits=5),gn,sh,alxSummStat)
end

function sampledAlpha(;d::Array,afs::Array,λdiv::Array,λpol::Array)

	#=pn = λpol[:,2]
	ps = λpol[:,1]=#

	## Outputs
	alphas, expDn, expDs = poissonFixation(observedValues=d,λds=λdiv[1],λdn=λdiv[2],λweak=λdiv[3],λstrong=λdiv[4],hri=hri)
	expPn, expPs         = poissonPolymorphism(observedValues=afs,λps=λpol[1],λpn=λpol[2])

	## Alpha from expected values. Used as summary statistics
	αS = @. round(1 - ((expDs/expDn) * (expPn/expPs)'),digits=5)

	return αS,alphas,expDn,expDs,expPn,expPs
end

function poissonFixation(;observedValues::Array, λds::Array, λdn::Array,λweak::Array,λstrong::Array,hri::Array)

    ds        = @. λds / (λds + λdn)
    dn        = @. λdn / (λds + λdn)
    dweak     = @. λweak / (λds + λdn)
    dstrong   = @. λstrong / (λds + λdn)
    dpositive = @. (λweak + λstrong) / (λds + λdn)
    dnohri    = @. ((λweak + λstrong)/hri) / (λds + λdn)

	sampledDs     = PoissonRandom.pois_rand.(ds .* observedValues)
	sampledDn     = PoissonRandom.pois_rand.(dn .* observedValues)
	sampledWeak   = PoissonRandom.pois_rand.(dweak .* observedValues)
	sampledStrong = PoissonRandom.pois_rand.(dstrong .* observedValues)
	sampledPositive = PoissonRandom.pois_rand.(dpositive .* observedValues)
	sampledNoHri = PoissonRandom.pois_rand.(dnohri .* observedValues)

	h = sampledPositive./sampledNoHri

	h = 1 .- ifelse.(h .> 1,1,h)

	alphas = @. [sampledWeak/sampledDn sampledStrong/sampledDn sampledPositive/sampledDn h]

	out = alphas,sampledDn, sampledDs
	return out
end
