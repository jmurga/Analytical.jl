################################
###    Summary statistics    ###
################################
"""
	poissonFixation(observedValues,λds, λdn)

Divergence sampling from Poisson distribution. The expected neutral and selected fixations are subset through their relative expected rates ([`Analytical.fixNeut`](@ref), [`Analytical.fixNegB`](@ref), [`Analytical.fixPosSim`](@ref)). Empirical values are used are used to simulate the locus *L* along a branch of time *T* from which the expected *Ds* and *Dn* raw count estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily given the success rate ``\\lambda`` in the distribution.

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
function poissonFixation(;observedValues::Array, λds::Float64, λdn::Float64)

	ds = λds / (λds + λdn)
	dn = λdn / (λds + λdn)
	observedValues = repeat(observedValues,1,10)
	# poissonS  = (ds .* observedValues) .|> Poisson
	# poissonD  = (dn .* observedValues) .|> Poisson
	# sampledDs = rand.(poissonS,1)
	# sampledDn = rand.(poissonD,1)

	sampledDs  = PoissonRandom.pois_rand.(ds .* observedValues)
	sampledDn  = PoissonRandom.pois_rand.(dn .* observedValues)

	out = sampledDn, sampledDs
	return out
end

"""
	poissonPolymorphism(observedValues,λps,λpn)

Polymorphism sampling from Poisson distributions. The total expected neutral and selected polimorphism are subset through the relative expected rates at the frequency spectrum ([`Analytical.fixNeut`](@ref), [`Analytical.DiscSFSNeutDown`](@ref),). Empirical sfs are used to simulate the locus *L* along a branch of time *T* from which the expected *Ps* and *Pn* raw count are estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily from the whole sfs given each frequency success rate ``\\lambda`` in the distribution.

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
function poissonPolymorphism(;observedValues::Array, λps::Array{Float64,1}, λpn::Array{Float64,1})

	λ1 = similar(λps);λ2 = similar(λpn)

	observedValues = repeat(observedValues,1,10)
	sampledPs = similar(observedValues)
	sampledPn = similar(observedValues)

	# Neutral λ
	λ1 .= @. λps / (λps + λpn)
	# λ1 = sum(λps) / (sum(λps) + sum(λpn))
	# Selected λ
	λ2 .= @. λpn / (λps + λpn)
	# λ2 = sum(λpn) / (sum(λps) + sum(λpn))

	# psPois(x,z=λ1) = rand.((z .* x ).|> Poisson,1)
	# pnPois(x,z=λ2) = rand.((z .* x ).|> Poisson,1)
	psPois(x,z=λ1) = PoissonRandom.pois_rand.((z .* x ))
	pnPois(x,z=λ2) = PoissonRandom.pois_rand.((z .* x ))

	sampledPs = psPois(observedValues)
	sampledPn = pnPois(observedValues)

	return (sampledPn, sampledPs)
	
end

"""
	sampledAlpha(observedValues,λds, λdn)

Ouput the expected values from the Poisson sampling process. Please check [`Analytical.poissonFixation`](@ref) and [`Analytical.poissonPolymorphism`](@ref) to understand the samplingn process. α(x) is estimated through the expected values of Dn, Ds, Pn and Ps.

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
function sampledAlpha(;param::parameters,d::Array,afs::Array,λdiv::Array{Float64,2},λpol::Array{Float64,2},bins::Int64=20)

	pn = λpol[:,2]
	ps = λpol[:,1]

	## Outputs
	expDn, expDs    = poissonFixation(observedValues=d,λds=λdiv[1],λdn=λdiv[2])
	expPn, expPs    = poissonPolymorphism(observedValues=afs,λps=ps,λpn=pn)

	if bins == (param.nn - 1)
		reducedExpPn = view(expPn,1:bins,:)
		reducedExpPs = view(expPs,1:bins,:)
	else
		reducedExpPn = view(reduceSfs(expPn,bins)',1:bins,:)
		reducedExpPs = view(reduceSfs(expPs,bins)',1:bins,:)
	end


	## Alpha from expected values. Used as summary statistics
	#ssAlpha = @. 1 - ((expDs/expDn) * (reducedExpPn./cumulativeExpPs))
	#ssAlpha = round.(ssAlpha,digits=5)

	αS = @. round(1 - ((expDs/expDn) * (expPn/expPs)),digits=5)

	return αS,expDn,expDs,expPn,expPs
end

"""
	alphaByFrequencies(gammaL,gammaH,pposL,param.pposH,data,nopos)

Analytical α(x) estimation. Solve α(x) from the expectation generally. We used the expected rates of divergence and polymorphism to approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region.

```math
\\mathbb{E}[\\alpha_{x}] =  1 - \\left(\\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{N}]}\\frac{\\mathbb{E}[P_{N}]}{\\mathbb{E}[P_{S}]}\\right)
```

# Arguments
 - `gammaL::Int64`: strength of weakly positive selection
 - `gammaH::Int64`: strength of strong positive selection
 - `pposL`::Float64: probability of weakly selected allele
 - `param.pposH`::Float64: probability of strong selected allele
 - `data::Array{Any,1}`: Array containing the total observed divergence, polymorphism and site frequency spectrum.
 - `nopos::String("pos","nopos","both")`: string to perform α(x) account or not for both positive selective alleles.

# Returns
 - `Array{Float64,1}` α(x).
"""
function analyticalAlpha(;param::parameters)

	##############################################################
						# Solve the model  #
	##############################################################
	B = param.B

	set_theta_f!(param)
	theta_f = param.theta_f
	# Solve the probabilities of fixations without background selection
	## First set non-bgs
	param.B = 0.999
	## Solve the mutation rate
	set_theta_f!(param)
	## Solve the probabilities
	setPpos!(param)
	# Return to the original values
	param.theta_f = theta_f
	param.B = B

	##############################################################
	# Accounting for positive alleles segregating due to linkage #
	##############################################################

	# Fixation
	fN     = param.B*fixNeut(param)
	fNeg   = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
	fPosL  = fixPosSim(param,param.gL,0.5*param.pposL)
	fPosH  = fixPosSim(param,param.gH,0.5*param.pposH)

	ds = fN
	dn = fNeg + fPosL + fPosH

	## Polymorphism
	neut = DiscSFSNeutDown(param,param.bn[param.B])

	selH = DiscSFSSelPosDown(param,param.gH,param.pposH,param.bn[param.B])
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL,param.bn[param.B])
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,param.bn[param.B])
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN))

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN

	# ps = @. neut / (sel+neut)
	# pn = @. sel / (sel+neut)

	## Outputs
	α = @. 1 - ((ds/dn) * (sel/neut))

	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	fN_nopos     = fN*(param.theta_mid_neutral/2.)*param.TE*param.NN
	fNeg_nopos   = fNeg*(param.theta_mid_neutral/2.)*param.TE*param.NN
	fPosL_nopos  = fPosL*(param.theta_mid_neutral/2.)*param.TE*param.NN
	fPosH_nopos  = fPosH*(param.theta_mid_neutral/2.)*param.TE*param.NN

	ds_nopos = fN_nopos
	dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

	## Polymorphism
	sel_nopos = selN
	ps_nopos = @. neut / (sel_nopos + neut)
	pn_nopos = @. sel_nopos / (sel_nopos + neut)

	α_nopos = 1 .- (ds_nopos/dn_nopos) .* (sel_nopos./neut)

	##########
	# Output #
	##########

	return (α,α_nopos)
	# return (α,α_nopos,[α[end] asymp1[1] c1[1] c2[1] c3[1]],[α_nopos[end] asymp2[1] c1[2] c2[2] c3[2]])
end

"""
	alphaByFrequencies(gammaL,gammaH,pposL,pposH,observedData,nopos)

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
function alphaByFrequencies(param::parameters,divergence::Array,sfs::Array,bins::Int64,cutoff::Float64,dac::Array{Float64,1})

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
	# neut = DiscSFSNeutDown(param,param.bn[param.B])
	neut = param.neut[param.B]

	selH = DiscSFSSelPosDown(param,param.gH,param.pposH,param.bn[param.B])
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL,param.bn[param.B])
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,param.bn[param.B])
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN))
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN
	
	## Outputs
	α = @. 1 - (ds/dn) * (sel/neut)
	# α = view(α,1:trunc(Int64,param.nn*cutoff),:)

	alxSummStat, expectedDn, expectedDs, expectedPn, expectedPs = sampledAlpha(param=param,d=divergence,afs=sfs,λdiv=hcat(ds,dn),λpol=hcat(neut,sel),bins=bins)

	# d=divergence;afs=sfs;λdiv=hcat(ds,dn);λpol=hcat(neut,sel);bins=bins
	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	fN_nopos       = fN*(param.theta_mid_neutral/2.)*param.TE*param.NN
	fNeg_nopos     = fNeg*(param.theta_mid_neutral/2.)*param.TE*param.NN
	fPosL_nopos    = fPosL*(param.theta_mid_neutral/2.)*param.TE*param.NN
	fPosH_nopos    = fPosH*(param.theta_mid_neutral/2.)*param.TE*param.NN

	ds_nopos       = fN_nopos
	dn_nopos       = fNeg_nopos + fPosL_nopos + fPosH_nopos
	dnS_nopos      = dn_nopos - fPosL_nopos

	## Polymorphism
	sel_nopos = selN

	## Outputs
	αW         = param.alLow/param.alTot
	α_nopos    = @. 1 - (ds_nopos/dn_nopos) * (sel_nopos/neut)[dac]
	# α_nopos    = @. 1 - (ds_nopos/dn_nopos) * (sel_nopos/neut)
	amk,ci,model    = asympFit(α_nopos)
	αW_nopos   = amk * αW
	αS_nopos   = amk * (1 - αW)
	# αW_nopos   = α_nopos * αW
	# αS_nopos   = α_nopos * (1 - αW)

	##########
	# Output #
	##########
	# Dn,Ds,Pn,Ps = expectedDn,expectedDs,sum(view(expectedPn,1,:),dims=2),sum(view(expectedPs,1,:),dims=2)

	# alphas = round.(hcat(param.alTot - param.alLow, param.alLow, param.alTot),digits=5)
	# alphas = round.(hcat(αW_nopos[trunc(Int64,param.nn*cutoff),:], αS_nopos[trunc(Int64,param.nn*cutoff),:], α_nopos[trunc(Int64,param.nn*cutoff),:]),digits=5)
	# alphas = round.(hcat(αW_nopos[end], αS_nopos[end], α_nopos[end]),digits=5)
	alphas = round.(hcat(αW_nopos, αS_nopos, amk),digits=5)
	alphas = repeat(alphas,outer=[10,1])

	# expectedValues = hcat(DataFrame(alphas),DataFrame(hcat(Dn,Ds,Pn,Ps)),DataFrame(permutedims(alxSummStat)),makeunique=true)
	expectedValues = hcat(alphas,permutedims(alxSummStat))

	return (α,α_nopos,expectedValues)
end

function writeSummaryStatistics(fileName::String,summStat)

	names = collect('a':'z')
	for i in 1:size(summStat,1)
		write(fileName * "_" * names[i] * ".tsv", DataFrame(summStat[i:i,:]), delim='\t', append=true)
	end
end
