################################
###    Summary statistics    ###
################################
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

function poissonPolymorphism2(;observedValues, λps, λpn)

    psPois(x,y=λps,z=λpn) = reduce(vcat,rand.((y./(y .+ z) .* x) .|> Poisson,1))
    pnPois(x,y=λps,z=λpn) = reduce(vcat,rand.((z./(y .+ z) .* x) .|> Poisson,1))

    sampledPs = observedValues .|> psPois # We can apply here any statistic measure
    sampledPn = observedValues .|> pnPois # We can apply here any statistic measure

    return (sum.(sampledPs), sum.(sampledPn))
end

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
		selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
		selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))

		sel = (selH+selL)+selN
		replace!(sel,NaN=>0)
		ps = sum(neut) ./ sum(sel+neut)
		pn = sum(sel) ./ sum(sel+neut)

		# Outputs
		expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
		expectedPs, expectedPn = poissonPolymorphism2(observedValues=[SFS],λps=ps,λpn=pn)


		sel = (selH+selL)+selN

		α = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)
		α = α[1:lastindex(α)-1]

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
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))
		sel = selN

		ps = neut ./ (sel.+neut)
		pn = sel ./ (sel.+neut)

		# Outputs
		expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
		expectedPs, expectedPn = poissonPolymorphism2(observedValues=[SFS],λps=ps,λpn=pn)

		α = 1 .- (fN/(fPosL + fPosH+  fNeg+0.0)) .* (sel./neut)
		α = α[1:lastindex(α)-1]

		expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,α[lastindex(α)])
		# expectedValues = hcat(expectedDs,expectedDn,expectedPs,expectedPn,fill(adap.B,size(expectedPn)[1]),fill(ret[lastindex(ret)],size(expectedPn)[1]))

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
		selH = cumulativeSfs(DiscSFSSelPosDown(gammaH,pposH))
		selL = cumulativeSfs(DiscSFSSelPosDown(gammaL,pposL))
		selN = cumulativeSfs(DiscSFSSelNegDown(pposH+pposL))

		sel = (selH+selL)+selN
		replace!(sel,NaN=>0)
		ps = sum(neut) ./ sum(sel+neut)
		pn = sum(sel) ./ sum(sel+neut)

		## Outputs
		expectedDs, expectedDn = poissonFixation(observedValues=D,λds=ds,λdn=dn)
		expectedPs, expectedPn = poissonPolymorphism2(observedValues=P,λps=ps,λpn=pn)

		α = 1 .- (fN/(fPosL + fPosH +  fNeg + 0.0)) .* (sel./neut)
		α = α[1:lastindex(α)-1]

		# Accounting only for neutral and deleterious alleles segregating
		## Fixation
		fN_nopos     = fN*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fNeg_nopos   = fNeg*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fPosL_nopos  = fPosL*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
		fPosH_nopos  = fPosH*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN

		ds_nopos = fN_nopos
		dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

		## Polymorphism
		sel_nopos = selN
		ps_nopos = sum(neut) ./ sum(sel_nopos+neut)
		pn_nopos = sum(sel_nopos) ./ sum(sel_nopos+neut)

		## Outputs
		expectedDs_nopos, expectedDn_nopos = poissonFixation(observedValues=D,λds=ds_nopos,λdn=dn_nopos)
		expectedPs_nopos, expectedPn_nopos = poissonPolymorphism2(observedValues=P,λps=ps_nopos,λpn=pn_nopos)

		α_nopos = 1 .- (fN_nopos/(fNeg_nopos + fPosH_nopos+  fNeg_nopos+0.0)) .* (sel_nopos./neut)
		α_nopos = α_nopos[1:lastindex(α_nopos)-1]

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
