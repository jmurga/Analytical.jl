################################
###    Summary statistics    ###
################################
"""
	poissonFixation(observedValues,λds, λdn)

Divergence sampling from Poisson distribution. The expected neutral and selected fixations are subset through their relative expected rates ([`Analytical.fixNeut`](@ref), [`Analytical.fixNegB`](@ref), [`Analytical.fixPosSim`](@ref)). Empirical values are used are used to simulate the locus *L* along a branch of time *T* from which the expected *Ds* and *Dn* raw count estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily given the success rate ``\\lambda`` in the distribution.

# Arguments
 - `observedValues::Array{Int64,1}`: Array containing the total observed divergence.
 - ` λds`: expected neutral fixations rate.
 - ` λdn`: expected selected fixations rate.
# Returns
 - `Array{Int64,1}` containing the expected count of neutral and selected fixations.

"""
function poissonFixation(;observedValues, λds, λdn)

	poissonS  = (λds/(λds + λdn) .* observedValues) .|> Poisson
	poissonD  = (λdn/(λds + λdn) .* observedValues) .|> Poisson

	sampledDs = rand.(poissonS,1)
	sampledDn = rand.(poissonD,1)

	return (reduce(hcat,sampledDs),reduce(hcat,sampledDn))
end
"""
	poissonPolymorphism(observedValues,λps,λpn)

Polymorphism sampling from Poisson distributions. The total expected neutral and selected polimorphism are subset through the relative expected rates at the frequency spectrum ([`Analytical.fixNeut`](@ref), [`Analytical.DiscSFSNeutDown`](@ref),). Empirical sfs are used to simulate the locus *L* along a branch of time *T* from which the expected *Ps* and *Pn* raw count are estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily from the whole sfs given each frequency success rate ``\\lambda`` in the distribution.

# Arguments
 - `observedValues::Array{Int64,1}`: Array containing the total observed divergence.
 - ` λps `: expected neutral site frequency spectrum rate.
 - ` λpn `: expected selected site frequency spectrum rate.
# Returns
 - `Array{Int64,1}` containing the expected total count of neutral and selected polymorphism.

"""
function poissonPolymorphism(;observedValues, λps, λpn)

	λ1 = similar(λps);λ2 = similar(λpn)
	sampledPs = similar(observedValues)
	sampledPn = similar(observedValues)
	
	# Neutral λ
	λ1 .= @. λps ./ (λps .+ λpn)
	# Selected λ
	λ2 .= @. λpn / (λps + λpn)

	psPois(x,z=λ1) = reduce.(vcat,rand.((z .* x ).|> Poisson,1))
	pnPois(x,z=λ2) = reduce.(vcat,rand.((z .* x ).|> Poisson,1))

	sampledPs = observedValues |> psPois
	sampledPn = observedValues |> pnPois

    return (sampledPs, sampledPn)
end


"""
	alphaByFrequencies(gammaL,gammaH,pposL,pposH,observedData,nopos)

Analytical α(x) estimation. Solve α(x) from the expectation generally. We used the expected rates of divergence and polymorphism to approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region.

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
 - `Array{Float64,1}` α(x).
"""
function analyticalAlpha(;gammaL::Int64,gammaH::Int64,pposL::Float64,pposH::Float64)

	##############################################################
	# Accounting for positive alleles segregating due to linkage #
	##############################################################

	# Fixation
	fN     = adap.B*fixNeut()
	fNeg   = adap.B*fixNegB(0.5*pposH+0.5*pposL)
	fPosL  = fixPosSim(gammaL,0.5*pposL)
	fPosH  = fixPosSim(gammaH,0.5*pposH)

	ds = fN
	dn = fNeg + fPosL + fPosH

	## Polymorphism
	neut = Array{Float64}(undef,adap.nn-1,1)
	selH = similar(neut);selL = similar(neut);selN = similar(neut);

	neut .= DiscSFSNeutDown()

	selH .= DiscSFSSelPosDown(gammaH,pposH)
	selL .= DiscSFSSelPosDown(gammaL,pposL)
	selN .= DiscSFSSelNegDown(pposH+pposL)
	splitColumns(matrix) = (view(matrix, :, i) for i in 1:size(matrix, 2))
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN))

	neut, selH, selL, selN = splitColumns(tmp)

	sel = (selH+selL)+selN
	sel = view(sel,1:lastindex(sel)-1,:)
	neut = view(neut,1:lastindex(neut)-1,:)

	if (isnan(sel[1]))
		sel[1]=0
	elseif(isnan(sel[end]))
		sel[end]=0
	end

	ps = similar(neut); pn = similar(neut)
	ps .= @. neut / (sel+neut)
	pn .= @. sel / (sel+neut)

	## Outputs
	α = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)

	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	fN_nopos     = fN*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
	fNeg_nopos   = fNeg*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
	fPosL_nopos  = fPosL*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
	fPosH_nopos  = fPosH*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN

	ds_nopos = fN_nopos
	dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

	## Polymorphism
	sel_nopos = view(selN,1:lastindex(selN)-1,:)
	ps_nopos  = similar(neut); pn_nopos = similar(neut)
	ps_nopos .= @. neut / (sel_nopos + neut)
	pn_nopos .= @. sel_nopos / (sel_nopos + neut)

	α_nopos = 1 .- (fN_nopos/(fPosL_nopos + fPosH_nopos +  fNeg_nopos + 0.0)) .* (sel_nopos./neut)

	##########
	# Output #
	##########
	return (α,α_nopos)
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
function alphaByFrequencies(;gammaL::Int64,gammaH::Int64,pposL::Float64,pposH::Float64,observedData)

	P   = observedData[1]
	sfs = observedData[2]
	D   = observedData[3]

	##############################################################
	# Accounting for positive alleles segregating due to linkage #
	##############################################################

	# Fixation
	fN     = adap.B*fixNeut()
	fNeg   = adap.B*fixNegB(0.5*pposH+0.5*pposL)
	fPosL  = fixPosSim(gammaL,0.5*pposL)
	fPosH  = fixPosSim(gammaH,0.5*pposH)

	ds = fN
	dn = fNeg + fPosL + fPosH

	## Polymorphism
	neut = Array{Float64}(undef,adap.nn-1,1)
	selH = similar(neut);selL = similar(neut);selN = similar(neut);

	neut .= DiscSFSNeutDown()

	selH .= DiscSFSSelPosDown(gammaH,pposH)
	selL .= DiscSFSSelPosDown(gammaL,pposL)
	selN .= DiscSFSSelNegDown(pposH+pposL)
	splitColumns(matrix) = (view(matrix, :, i) for i in 1:size(matrix, 2))
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN))

	neut, selH, selL, selN = splitColumns(tmp)

	sel = (selH+selL)+selN
	sel = view(sel,1:lastindex(sel)-1,:)
	neut = view(neut,1:lastindex(neut)-1,:)

	if (isnan(sel[1]))
		sel[1]=0
	elseif(isnan(sel[end]))
		sel[end]=0
	end

	ps = similar(neut); pn = similar(neut)
	ps .= @. neut / (sel+neut)
	pn .= @. sel / (sel+neut)

	## Outputs
	expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
	expectedPs, expectedPn = poissonPolymorphism(observedValues=sfs,λps=ps,λpn=pn)

	cumulativePs = view(cumulativeSfs(expectedPs),1:size(expectedPn,1),:)
	cumulativePn = view(cumulativeSfs(expectedPn),1:size(expectedPn,1),:)

	# α = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)
	α = view(1 .- (((expectedDs)./(expectedDn)) .* (cumulativePn./cumulativePs)),1:convert(Int64,ceil(adap.nn*0.9)),:)

	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	fN_nopos     = fN*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
	fNeg_nopos   = fNeg*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
	fPosL_nopos  = fPosL*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
	fPosH_nopos  = fPosH*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN

	ds_nopos = fN_nopos
	dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

	## Polymorphism
	sel_nopos = view(selN,1:lastindex(selN)-1,:)
	ps_nopos  = similar(neut); pn_nopos = similar(neut)
	ps_nopos .= @. neut / (sel_nopos + neut)
	pn_nopos .= @. sel_nopos / (sel_nopos + neut)

	## Outputs
	expectedDs_nopos, expectedDn_nopos = poissonFixation(observedValues=D,λds=ds_nopos,λdn=dn_nopos)
	expectedPs_nopos, expectedPn_nopos = poissonPolymorphism(observedValues=sfs,λps=ps_nopos,λpn=pn_nopos)

	cumulativePs_nopos = view(cumulativeSfs(expectedPs_nopos),1:adap.nn-1,:)
	cumulativePn_nopos = view(cumulativeSfs(expectedPn_nopos),1:adap.nn-1,:)

	# α_nopos = 1 .- (fN_nopos/(fPosL_nopos + fPosH_nopos +  fNeg_nopos + 0.0)) .* (sel_nopos./neut)
	α_nopos = view(1 .- ((expectedDs_nopos ./ expectedDn_nopos) .* (cumulativePn_nopos ./ cumulativePs_nopos)),1:convert(Int64,ceil(adap.nn*0.9)),:)

	##########
	# Output #
	##########
	
	summarySfs = reduceSfs(expectedPn + expectedPn,20)

	# Handling error to return any array size
	Dn,Ds,Pn,Ps = try 
		permutedims(expectedDs),permutedims(expectedDn),permutedims(sum(expectedPs,dims=1)),permutedims(sum(expectedPn,dims=1))
	catch err
		expectedDs,expectedDn,sum(expectedPs,dims=1),sum(expectedPn,dims=1)
	end

	expectedValues = hcat(Dn,Ds,Pn,Ps,summarySfs,view(α,size(α,1),:), view(α_nopos,size(α_nopos,1),:) .- view(α,size(α,1),:), view(α_nopos,size(α_nopos,1),:))

	return (α,α_nopos,expectedValues)
end

function summaryStatistics(fileName::String,summStats::Array{Float64})

	write(fileName, DataFrame(summStats), delim='\t', append=true)

end
