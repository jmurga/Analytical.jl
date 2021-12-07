################################
###    Summary statistics    ###
################################
"""
	poisson_fixation(empirical_values,λds, λdn)

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
function poisson_fixation(;empirical_values::Array, λds::Array, λdn::Array,λweak::Array,λstrong::Array)

    ds            = @. λds / (λds + λdn) * empirical_values
    dn            = @. λdn / (λds + λdn) * empirical_values
    dn_weak         = @. λweak / (λds + λdn) * empirical_values
    dn_strong       = @. λstrong / (λds + λdn) * empirical_values

    sampled_Ds     = pois_rand.(ds)
    sampled_Dn     = pois_rand.(dn)
    sampled_weak   = pois_rand.(dn_weak)
    sampled_strong = pois_rand.(dn_strong)

    alphas        = @. [sampled_weak/sampled_Dn sampled_strong/sampled_Dn (sampled_weak+sampled_strong)/sampled_Dn]

	out = alphas,sampled_Dn, sampled_Ds
	return out
end

"""
	poisson_polymorphism(empirical_values,λps,λpn)

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
function poisson_polymorphism(;empirical_values::Array, λps::Array, λpn::Array,pol::Array)

	# Neutral λ;
	λ1 = @. λps / (λps + λpn) * empirical_values'
	# Selected λ;
	λ2 = @. λpn / (λps + λpn) * empirical_values'
	
	# Neutral λ;
	#=λ1 = @. λps / pol[:,1] * observedValues[:,2]'
	# Selected λ;
	λ2 = @. λpn / pol[:,2] * observedValues[:,1]'=#

	# Replace negative and 0 rates values. Strange behaviour depending on θ and Γ values
	# Relative rates output NaN values due to 0 divisons.

	replace!(λ1,NaN=>1)	
	replace!(λ2,NaN=>1)

	sampled_Ps = pois_rand.(λ1)
	sampled_Pn = pois_rand.(λ2)

	return (sampled_Pn, sampled_Ps)
end

"""
	sampledAlpha(observedValues,λds, λdn)

Ouput the expected values from the Poisson sampling process. Please check [`poisson_fixation`](@ref) and [`poisson_polymorphism`](@ref) to understand the samplingn process. α(x) is estimated through the expected values of Dn, Ds, Pn and Ps.

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

"""

	Ouput the expected values from the Poisson sampling process. Please check [`poisson_fixation`](@ref) and [`poisson_polymorphism`](@ref) to understand the samplingn process. α(x) is estimated through the expected values of Dn, Ds, Pn and Ps.

	# Arguments
	 - `afs::Vector`: Array containing SFS subset
	 - `divergence::Vector`: Array containing divergence counts
	 - `observed::Matrix`: Array containing SFS
	 - `m::Array`: Array containing modeled parameters
	 - `nt::Array`: Array containing neutral polymorphic rates at frequencies dac
	 - `sl::Array`: Array containing selective polymorphic rates at frequencies dac
	 - `d::Array`: Array containing fixation rates
	 - `p::Array`: Array containing total polymorphic rates
	 - `dac::Array`: Array containing Derived Allel Counts to subset
	# Returns
	Expected α_x values
	 - `Array{Int64,2}` true α, γ, shape parameter and α_x values.

	sampling_from_rates(afs,divergence,observed,m,nt,d,p,dac)
"""
function sampling_from_rates(afs::Vector,divergence::Vector,observed::Matrix,m::Array,nt::Array,sl::Array,d::Array,p::Array,dac::Array)

	ds             = d[:,1]
	dn             = d[:,2]
	dn_weak          = d[:,3]
	dn_strong        = d[:,4]
	gn             = abs.(m[:,4])
	sh             = round.(m[:,end-1],digits=5)

	#alxSummStat, alphasDiv, expectedDn, expectedDs, expectedPn, expectedPs = sampledAlpha(divergence=divergence,afs=s,observed=observed,λdiv=[ds,dn,dn_weak,dn_strong],λpol=[permutedims(nt),permutedims(sl)],dac=dnac)

	## Outputs
	alphas, expDn, expDs = poisson_fixation(empirical_values=divergence,λds=ds,λdn=dn,λweak=dn_weak,λstrong=dn_strong)
	expPn, expPs         = poisson_polymorphism(empirical_values=afs,λps=nt,λpn=sl,pol=p)
	
	#=θ_pn = @. observed[1,1] / (expPn[1,:]);
	θ_ps = @. observed[1,2] / (expPs[1,:]);

	r_pn = ones(size(expPn));
	r_ps = ones(size(expPs));
	for j = 2:size(dac,1)
		r_pn[j,:] .= @. (observed[j,1]) / (expPn[j,:]);
		r_pn[j,:] .= @. r_pn[j,:] / θ_pn;
		
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
	end=#

	## Alpha from expected values. Used as summary statistics
	α_x = @. round(1 - ((expDs/expDn) * (expPn/expPs)),digits=5)
	#=α_x = @. round(1 - ((expDs/expDn) * (expPn_rj/expPs_rj)'),digits=5)=#


	expected_values = hcat(round.(alphas,digits=5),gn,sh,α_x)
	
	return(expected_values)
end

"""
	summary_statistics(param::parameters,rates::JLD2.JLDFile,analysis_folder::String,summstat_size::Int64,replicas::Int64,bootstrap::Bool)

Estimate summary statistics using observed data and analytical rates. *analysis_folder* will check for the SFS and divergence file and will be used to output summary statistics

# Arguments
 - `param::parameters` : Mutable structure containing the models
 - `rates::JLD2.JLDFile` : HDF5 containing solved models, fixation and polymorphic rates
 - `analysis_folder::String` : Folder containing the SFS and divergence files. It will be used to output the observed data and summary estatistics.
 - `summstat_size::Int64` : Number of summary statistics
 - `replicas::Int64` : Number of bootstrap replicas
 - `bootstrap::Bool` : Boolean to perform or not bootstrapping
# Output
 - Obserded data and summary statistics to ABC inference
"""
function summary_statistics(;param::parameters,h5_file::String,analysis_folder::String,summstat_size::Int64,replicas::Int64=1,bootstrap::Bool=false)

	#Opening files
    sfs_files        = filter(x -> occursin("sfs",x), readdir(analysis_folder,join=true));
    divergence_files = filter(x -> occursin("div",x), readdir(analysis_folder,join=true));

	sfs,divergence,α,observed = open_sfs_div(sfs_files,divergence_files,param.dac,replicas,bootstrap);

	#Open rates
	h5    = jldopen(h5_file)
	tmp   = h5[string(param.N) * "/" * string(param.n)]

	#Subset index
	idx    = sample(1:size(tmp["models"],1),summstat_size,replace=false)
	models = Array(view(tmp["models"],idx,:));
	dsdn   = Array(view(tmp["dsdn"],idx,:));
	pol    = Array(view(tmp["pol"],idx,:));

	# Filtering polymorphic rate by dac
	n    = hcat(map(x -> view(tmp["neut"][x],:),param.dac)...);
	s    = hcat(map(x -> view(tmp["sel"][x],:),param.dac)...);
	neut = Array(view(n,idx,:));
	sel  = Array(view(s,idx,:));

	#Making summaries
	f(afs,divergence,observed,m=models,nt=neut,sl=sel,d=dsdn,p=pol,dac=param.dac) = sampling_from_rates(afs,divergence,observed,m,nt,sl,d,p,dac);

	#expectedValues = f.(sfs,divergence,observed);
	expectedValues = progress_pmap(f,sfs,divergence,observed;progress= Progress(size(sfs,1),desc="Sampling P and D "));

	w(x,name) = CSV.write(name,DataFrame(x,:auto),delim='\t',header=false);

	# Controling outlier cases
	fltInf(e) = replace!(e, -Inf=>NaN)
	expectedValues = fltInf.(expectedValues)
	fltNan(e) = e[vec(.!any(isnan.(e),dims=2)),:]
	expectedValues = fltInf.(expectedValues)
	expectedValues = fltNan.(expectedValues)
	#=expectedValues = map(x -> x[(x[:,3] .> 0) .& (x[:,3] .<1 ),:],expectedValues)=#
	expectedValues = map(x -> x[(x[:,3] .< 1),:],expectedValues)
	
	# Writting ABCreg input
	#=@inbounds @sync for i in eachindex(α)
		Base.Threads.@spawn begin
			w(α,analysis_folder * "/alphas_" * string(i) * ".txt")
			w(expectedValues,  analysis_folder * "/summstat_" * string(i) * ".txt")
		end
	end=#
	progress_map(w,α, analysis_folder * "/alphas_" .* string.(1:size(sfs,1)) .* ".txt";progress= Progress(size(sfs_files,1),desc="Writting α "));
	progress_pmap(w,expectedValues,  analysis_folder * "/summstat_" .* string.(1:size(sfs,1)) .* ".txt";progress= Progress(size(sfs,1),desc="Writting summaries "));

	return(expectedValues)
end
