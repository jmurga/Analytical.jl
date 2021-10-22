################################
###    Summary statistics    ###
################################
"""
	poissonFixation(observedValues,λds, λdn)

Divergence sampling from Poisson distribution. The expected neutral and selected fixations are subset through their relative expected rates ([`fixNeut`](@ref), [`fixNegB`](@ref), [`fixPosSim`](@ref)). Empirical values are used are used to simulate the locus *L* along a branch of time *T* from which the expected *Ds* and *Dn* raw count estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily given the success rate ``\\lambda`` in the distribution.

```math
\\mathbb{E}[D_N] = X \\in Poisson\\left(\\lambda = D \\times \\left[\\frac{\\mathbb{E}[D_+] + \\mathbb{E}[D_-]}{\\mathbb{E}[D_+] + \\mathbb{E}[D_-] + \\mathbb{E}[D_0]}\\right]\\right)
```
```math
\\mathbb{E}[D_S] = X \\in Poisson\\left(\\lambda = D \\times \\left[\\frac{\\mathbb{E}[D_0]}{\\mathbb{E}[D_+] + \\mathbb{E}[D_-] + \\mathbb{E}[D_0]}\\right]\\right)
```
# Arguments
 - `observedValues::Array`: Array containing the total observed divergence.
 - ` λds::Float64`: expected neutral fixations rate.
 - ` λdn::Float64`: expected selected fixations rate.
# Returns
 - `Array{Int64,1}` containing the expected count of neutral and selected fixations.

"""
function poissonFixation(;observedValues::Array, λds::Array, λdn::Array,λweak::Array,λstrong::Array)

	ds = @. λds / (λds + λdn) * observedValues
	dn = @. λdn / (λds + λdn) * observedValues
	dweak = @. λweak / (λds + λdn) * observedValues
	dstrong = @. λstrong / (λds + λdn) * observedValues

	sampledDs     = pois_rand.(ds)
	sampledDn     = pois_rand.(dn)
	sampledWeak   = pois_rand.(dweak)
	sampledStrong = pois_rand.(dstrong)

	alphas = @. [sampledWeak/sampledDn sampledStrong/sampledDn (sampledWeak+sampledStrong)/sampledDn]

	out = alphas,sampledDn, sampledDs
	return out
end

"""
	poissonPolymorphism(observedValues,λps,λpn)

Polymorphism sampling from Poisson distributions. The total expected neutral and selected polimorphism are subset through the relative expected rates at the frequency spectrum ([`fixNeut`](@ref), [`DiscSFSNeutDown`](@ref),). Empirical sfs are used to simulate the locus *L* along a branch of time *T* from which the expected *Ps* and *Pn* raw count are estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily from the whole sfs given each frequency success rate ``\\lambda`` in the distribution.

The success rate managing the Poisson distribution by the observed count each frequency.  We considered both sampling variance and process variance is affecting the number of variable alleles we sample from SFS. This variance arises from the random mutation-fixation process along the branch. To incorporate this variance we do one sample per frequency-bin and use the total sampled variation and the SFS for the summary statistics.

```math
\\mathbb{E}[P_N] = \\sum_{x=0}^{x=1} X \\in Poisson\\left(\\lambda = SFS_{(x)} \\times \\left[\\frac{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}]}{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}] + \\mathbb{E}[P_{0(x)}]}\\right]\\right)
```

```math
\\mathbb{E}[P_S] = \\sum_{x=0}^{x=1} X \\in Poisson\\left(\\lambda = SFS_{(x)} \\times \\left[\\frac{\\mathbb{E}[P_{0(x)}]}{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}] + \\mathbb{E}[P_{0(x)}]}\\right]\\right)
```

# Arguments
 - `observedValues::Array{Int64,1}`: Array containing the total observed divergence.
 - ` λps::Array{Float64,1} `: expected neutral site frequency spectrum rate.
 - ` λpn::Array{Float64,1} `: expected selected site frequency spectrum rate.
# Returns
 - `Array{Int64,2}` containing the expected total count of neutral and selected polymorphism.

"""
function poissonPolymorphism(;observedValues::Array, λps::Array, λpn::Array)

	# Neutral λ;
	λ1 = @. λps / (λps + λpn) * observedValues
	# Selected λ;
	λ2 = @. λpn / (λps + λpn) * observedValues

	# Replace negative and 0 rates values. Strange behaviour depending on θ and Γ values
	# Relative rates output NaN values due to 0 divisons.

	replace!(λ1,NaN=>1)	
	replace!(λ2,NaN=>1)

	sampledPs = pois_rand.(λ1)
	sampledPn = pois_rand.(λ2)

	return (sampledPn, sampledPs)
end

"""
	sampledAlpha(observedValues,λds, λdn)

Ouput the expected values from the Poisson sampling process. Please check [`poissonFixation`](@ref) and [`poissonPolymorphism`](@ref) to understand the samplingn process. α(x) is estimated through the expected values of Dn, Ds, Pn and Ps.

# Arguments
 - `param::parameters`: Array containing the total observed divergence.
 - `d::Array`: observed divergence.
 - `afs::Array`: observed polymorphism.
 - ` λdiv::Array{Float64,2}`: expected fixations rate.
 - ` λdiv::Array{Float64,2}`: expected site frequency spectrum rates.
# Returns
αS,expDn,expDs,expPn,expPs,ssAlpha
 - `Array{Int64,2}` containing α(x) values.
 - `Array{Int64,1}` expected non-synonymous divergence.
 - `Array{Int64,1}` expected synonymous divergence.
 - `Array{Int64,1}` expected non-synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,2}` containing α(x) binned values.

"""
function sampledAlpha(;d::Array,afs::Array,observed::Matrix,λdiv::Array,λpol::Array,dac::Array)

	## Outputs
	alphas, expDn, expDs = poissonFixation(observedValues=d,λds=λdiv[1],λdn=λdiv[2],λweak=λdiv[3],λstrong=λdiv[4])
	expPn, expPs         = poissonPolymorphism(observedValues=afs,λps=λpol[1],λpn=λpol[2])
	
	θ_pn = @. observed[1,1] / (expPn[1,:]);
	θ_ps = @. observed[1,2] / (expPs[1,:]);

	r_pn = ones(size(expPn));
	r_ps = ones(size(expPs));
	for j = 2:size(dac,1)
		r_pn[j,:] .= @. (observed[j,1]) / (expPn[j,:]);
		r_pn[j,:] = @. r_pn[j,:] / θ_pn;
		
		r_ps[j,:] .= @. (observed[j,2]) / (expPs[j,:]);
		r_ps[j,:] .= @. r_ps[j,:] / θ_ps;


	end

	expPn_rj = zeros(size(expPn))
	expPs_rj = zeros(size(expPs))

	expPn_rj[1,:] = @. expPn[1,:] * θ_pn;
	expPs_rj[1,:] = @. expPs[1,:] * θ_ps;
	for j = 2:size(dac,1)
		expPn_rj[j,:] = @. expPn[j,:] * θ_pn * r_pn[j,:];
		expPs_rj[j,:] = @. expPs[j,:] * θ_ps * r_ps[j,:];
	end

	## Alpha from expected values. Used as summary statistics
	αS = @. round(1 - ((expDs/expDn) * (expPn/expPs)'),digits=5)
	αS_v2 = @. round(1 - ((expDs/expDn) * (expPn_rj/expPs_rj)'),digits=5)

	return αS,alphas,expDn,expDs,expPn,expPs
end

"""
	samplingFromRates(gammaL,gammaH,pposL,pposH,observedData,nopos)
"""
function samplingFromRates(m::Array,s::Array,d::Array,nt::Array,sl::Array,x::Array,observed::Matrix,dac::Array)

	ds             = x[:,1]
	dn             = x[:,2]
	dweak          = x[:,3]
	dstrong        = x[:,4]
	gn             = abs.(m[:,4])
	sh             = round.(m[:,end-1],digits=5)

	alxSummStat, alphasDiv, expectedDn, expectedDs, expectedPn, expectedPs = sampledAlpha(d=d,afs=s,observed=observed,λdiv=[ds,dn,dweak,dstrong],λpol=[permutedims(nt),permutedims(sl)],dac=dac)
	expectedValues = hcat(round.(alphasDiv,digits=5),gn,sh,alxSummStat)
	
	return(expectedValues)
end

"""
	summaryStatsFromRates(param::parameters,rates::JLD2.JLDFile,analysisFolder::String,summstatSize::Int64,replicas::Int64,bootstrap::Bool)

Estimate summary statistics using observed data and analytical rates. *analysisFolder* will check for the SFS and divergence file and will be used to output summary statistics

# Arguments
 - `param::parameters` : Mutable structure containing the models
 - `rates::JLD2.JLDFile` : HDF5 containing solved models, fixation and polymorphic rates
 - `analysisFolder::String` : Folder containing the SFS and divergence files. It will be used to output the observed data and summary estatistics.
 - `summstatSize::Int64` : Number of summary statistics
 - `replicas::Int64` : Number of bootstrap replicas
 - `bootstrap::Bool` : Boolean to perform or not bootstrapping
# Output
 - Obserded data and summary statistics to ABC inference
"""
function summaryStatsFromRates(;param::parameters,h5file::JLD2.JLDFile,analysisFolder::String,summstatSize::Int64,replicas::Int64,bootstrap::Bool)

	#Opening files
	sFile   = filter(x -> occursin("sfs",x), readdir(analysisFolder,join=true));
	dFile   = filter(x -> occursin("div",x), readdir(analysisFolder,join=true));

	sfs,d,α,observed = openSfsDiv(sFile,dFile,param.dac,replicas,bootstrap);

	#Open rates
	tmp    = h5file[string(param.N) * "/" * string(param.n)]

	#Subset index
	idx    = sample(1:size(tmp["models"],1),summstatSize,replace=false)
	models = Array(view(tmp["models"],idx,:));
	dsdn   = Array(view(tmp["dsdn"],idx,:));

	# Filtering polymorphic rate by dac
	n    = hcat(map(x -> view(tmp["neut"][x],:),param.dac)...);
	s    = hcat(map(x -> view(tmp["sel"][x],:),param.dac)...);
	neut = Array(view(n,idx,:));
	sel  = Array(view(s,idx,:));

	#Making summaries
	f(x,y,o,m=models,nt=neut,sl=sel,d=dsdn,dac=param.dac) = samplingFromRates(m,x,y,nt,sl,d,o,dac);
	expectedValues = f.(sfs,d,observed);

	w(x,name) = CSV.write(name,DataFrame(x,:auto),delim='\t',header=false);

	# Controling outlier cases
	fltInf(e) = replace!(e, -Inf=>NaN)
	expectedValues = fltInf.(expectedValues)
	fltNan(e) = e[vec(.!any(isnan.(e),dims=2)),:]
	expectedValues = fltInf.(expectedValues)
	expectedValues = fltNan.(expectedValues)
	expectedValues = map(x -> x[(x[:,3] .> 0 ) .& (x[:,3] .<1 ),:],expectedValues)
	
	# Writting ABCreg input
	#=w(vcat(α...), analysisFolder * "/alphas.txt");
	w(expectedValues, analysisFolder * "/summstat.txt");=#
	progress_map(w,α, analysisFolder * "/alphas_" .* string.(1:size(sFile,1)) .* ".txt";progress= Progress(size(sFile,1),desc="Writting α "));
	progress_pmap(w,expectedValues,  analysisFolder * "/summstat_" .* string.(1:size(sFile,1)) .* ".txt";progress= Progress(size(sFile,1),desc="Writting summaries "));

	return(expectedValues)
end
