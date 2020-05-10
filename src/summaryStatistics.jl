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

	return(reduce(vcat,sampledDs),reduce(vcat,sampledDn))
end

function poissonPolymorphism(;observedValues, λps, λpn)

    psPois(x,y=λps,z=λpn) = reduce(vcat,rand.((y./(y .+ z) .* x) .|> Poisson,1))
    pnPois(x,y=λps,z=λpn) = reduce(vcat,rand.((z./(y .+ z) .* x) .|> Poisson,1))

    sampledPs = observedValues .|> psPois # We can apply here any statistic measure
    sampledPn = observedValues .|> pnPois # We can apply here any statistic measure

    return sampledPs,sampledPn
end

"""
	poissonPolymorphism(observedValues,λps,λpn)

Polymorphism sampling from Poisson distributions. The total expected neutral and selected polimorphism are subset through the relative expected rates at the frequency spectrum ([`Analytical.fixNeut`](@ref), [`Analytical.DiscSFSNeutDown`](@ref),). Empirical SFS are used to simulate the locus *L* along a branch of time *T* from which the expected *Ps* and *Pn* raw count are estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily from the whole SFS given each frequency success rate ``\\lambda`` in the distribution.

# Arguments
 - `observedValues::Array{Int64,1}`: Array containing the total observed divergence.
 - ` λps `: expected neutral site frequency spectrum rate.
 - ` λpn `: expected selected site frequency spectrum rate.
# Returns
 - `Array{Int64,1}` containing the expected total count of neutral and selected polymorphism.

"""
function poissonPolymorphism2(;observedValues, λps, λpn)

    psPois(x,y=λps,z=λpn) = reduce(vcat,rand.((y./(y .+ z) .* x) .|> Poisson,1))
    pnPois(x,y=λps,z=λpn) = reduce(vcat,rand.((z./(y .+ z) .* x) .|> Poisson,1))

    sampledPs = observedValues .|> psPois # We can apply here any statistic measure
    sampledPn = observedValues .|> pnPois # We can apply here any statistic measure

    return (sum.(sampledPs), sum.(sampledPn))
end

"""
	alphaByFrequencies(gammaL,gammaH,pposL,pposH,observedData,nopos)

Analytical α(x) estimation. We used the expected rates of divergence and polymorphism to approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region. Solve α(x) from the expectation rates:

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
function alphaByFrequencies(gammaL::Int64,gammaH::Int64,pposL::Float64,pposH::Float64,observedData::Array,nopos::String)

	P   = observedData[1,:][lastindex(observedData[1,:])]
	SFS = observedData[2,:][lastindex(observedData[2,:])]
	D = observedData[3,:][lastindex(observedData[3,:])]

	if nopos == "pos"
		# Fixation
		fN    = adap.B*fixNeut()
		fNeg  = adap.B*fixNegB(0.5*pposH+0.5*pposL)
		fPosL = fixPosSim(gammaL,0.5*pposL)
		fPosH = fixPosSim(gammaH,0.5*pposH)

		ds = fN
		dn = fNeg + fPosL + fPosH

		# Polymorphism
		neut = cumulativeSfs(DiscSFSNeutDown())
		neut = neut[1:lastindex(neut)-1]

		selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
		selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))

		sel = (selH+selL)+selN
		sel = sel[1:lastindex(sel)-1]

		ps = neut ./ sum(sel.+neut)
		pn = sel ./ sum(sel.+neut)
		
		# Outputs
		expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
		expectedPs, expectedPn = poissonPolymorphism2(observedValues=[SFS],λps=ps,λpn=pn)

		α = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)

		expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,α[lastindex(α)])

		return (α,expectedValues)
	elseif nopos == "nopos"

		# Fixation
		fN    = adap.B*fixNeut()*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN
		fNeg  = adap.B*fixNegB(0.5*pposH+0.5*pposL)*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN
		fPosL = fixPosSim(gammaL,0.5*pposL)*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN
		fPosH = fixPosSim(gammaH,0.5*pposH)*(adap.theta_mid_neutral/2.0)*adap.TE*adap.NN

		ds = fN
		dn = fNeg + fPosL + fPosH

		# Polymorphism
		neut = cumulativeSfs(DiscSFSNeutDown())
		neut = neut[1:lastindex(neut)-1]

		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))
		sel = selN[1:lastindex(selN)-1]
		
		ps = neut ./ (sel.+neut)
		pn = sel ./ (sel.+neut)
		

		# Outputs
		expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
		expectedPs, expectedPn = poissonPolymorphism2(observedValues=[SFS],λps=ps,λpn=pn)

		α = 1 .- (fN/(fPosL + fPosH+  fNeg+0.0)) .* (sel./neut)

		expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,α[lastindex(α)])
		return (α,expectedValues)
	else

		## Accounting for positive alleles segregating due to linkage
		# Fixation
		fN     = adap.B*fixNeut()
		fNeg   = adap.B*fixNegB(0.5*pposH+0.5*pposL)
		fPosL  = fixPosSim(gammaL,0.5*pposL)
		fPosH  = fixPosSim(gammaH,0.5*pposH)

		ds = fN
		dn = fNeg + fPosL + fPosH

		## Polymorphism
		neut = cumulativeSfs(DiscSFSNeutDown())
		neut = neut[1:lastindex(neut)-1]
		
		selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
		selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))

		sel = (selH+selL)+selN
		sel = sel[1:lastindex(sel)-1]
		ps = neut ./ sum(sel.+neut)
		pn = sel ./ sum(sel.+neut)

		## Outputs
		expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
		expectedPs, expectedPn = poissonPolymorphism2(observedValues=P,λps=ps,λpn=pn)

		α = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)

		# Accounting only for neutral and deleterious alleles segregating
		## Fixation
		fN_nopos     = fN*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fNeg_nopos   = fNeg*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fPosL_nopos  = fPosL*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fPosH_nopos  = fPosH*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN

		ds_nopos = fN_nopos
		dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

		## Polymorphism
		sel_nopos = selN[1:lastindex(selN)-1]
		ps_nopos = neut ./ (sel_nopos.+neut)
		pn_nopos = sel_nopos ./ (sel_nopos.+neut)

		## Outputs
		expectedDs_nopos, expectedDn_nopos = poissonFixation(observedValues=D,λds=ds_nopos,λdn=dn_nopos)
		expectedPs_nopos, expectedPn_nopos = poissonPolymorphism2(observedValues=[SFS],λps=ps_nopos,λpn=pn_nopos)

		α_nopos = 1 .- (fN_nopos/(fPosL_nopos + fPosH_nopos +  fNeg_nopos + 0.0)) .* (sel_nopos./neut)

		expectedValues = hcat(expectedDs_nopos,expectedDn_nopos,expectedPs_nopos,expectedPn_nopos,α[lastindex(α)],α_nopos[lastindex(α_nopos)])

		return (α,α_nopos,expectedValues)
	end
end

function summaryStatistics(fileName,alpha,expectedValues)


	h5open(fileName, "cw") do file
		tmp = string(rand(Int64))
	    write(file, "tmp"*tmp*"/alpha", alpha)
	    write(file, "tmp"*tmp*"/expectedValues", expectedValues)
	end
end
