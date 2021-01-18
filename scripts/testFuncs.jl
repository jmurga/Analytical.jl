# Open empirical data
# param = parameters(N=1000,n=50,B=0.999,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=0.4,alLow=0.2);param.nn=101 ;binomOp(param)
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.2_0.999/sfs.tsv")))
# sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/tennesen/tennesen_0.4_0.1_0.999/sfs.tsv")))
sfs = sfs[:,2:end]

sCumu = convert.(Int64,Analytical.cumulativeSfs(sfs))

sfsPos   = sCumu[:,1] + sCumu[:,2]
sfsNopos = sCumu[:,4] + sCumu[:,2]

div = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.2_0.999/div.tsv")))
d = [convert(Int64,sum(divergence[1:2]))]

rSfs     = Analytical.reduceSfs(sfs,100)'
alpha    = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'
rSfsCumu    = Analytical.reduceSfs(sCumu,100)'
alphaCumu   = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)'


param = Analytical.parameters(N=1000,n=661,B=0.999,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=0.4,alLow=0.3,Lf=2*10^5)
Analytical.binomOp!(param)

@time j = param.B;
@time Analytical.set_theta_f!(param);
@time theta_f = param.theta_f;
@time param.B = 0.999;
@time Analytical.set_theta_f!(param);
@time Analytical.setPpos!(param);
@time param.theta_f = theta_f;
@time param.B = j;

x3 ,y3 = Analytical.analyticalAlpha(param=param)

x1,y1,z1 = alphaByFrequencies(param,d,s,b,c)
x2,y2,z2 = alphaByFrequencies(param,d,sum(sCumu,dims=2),b,c)

Plots.theme(:vibrant)
p = hcat(x,y,reverse(alphaCumu'),convert(Array,z2[:,4:end])')
lb=["Analytical α" "Analytical α_nopos" "Simulated α" "Sampled α"]
fn = plot(p,legend = :outertopright,label=lb,layout=4)
fn1 = plot(α,label="Analytical α",legend = :outertopright)
fn2 = plot(alxSummStat,label="Sampled α",legend = :outertopright,linecolor="#d43625")
l = @layout [a; b; c]
plot(fn,fn1,fn2,layout=l)

function summStats(param::parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
    # @threads

        @showprogress for i in 1:iter
        # for i in 1:iter


                fac       = rand(-2:0.05:2)
                afac      = 0.184*(2^fac)
                bfac      = 0.000402*(2^fac)

                # alTot     = rand(collect(0.01:0.01:0.4))
                alTot     = rand(collect(0.05:0.05:0.4))
                # alLow     = round(rand(collect((alTot/10):(alTot/10):alTot)),digits=5)
                lfac      = rand(collect(0.1:0.1:0.9))
                alLow     = round(alTot * lfac,digits=5)
        # println((thread=Threads.threadid(), iteration=i))

                bgsIter(param,afac,bfac,alTot,alLow,div,sfs,output,b,c)
        end
end

function bgsIter(param::parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)

        for j in param.bRange
                # j = 0.999
                param.al = afac; param.be = bfac;
                param.alLow = alLow; param.alTot = alTot; param.B = j

                set_theta_f!(param)
                theta_f = param.theta_f
                param.B = 0.999
                set_theta_f!(param)
                setPpos!(param)
                param.theta_f = theta_f
                param.B = j
                # x,y = analyticalAlpha(param=param)
                x,y,z = alphaByFrequencies(param,div,sfs,b,c)
                # println(z)

                summaryStatistics(output, z)

        end
end

# pol,sfs,div = parseSfs(param=param,data=files[1],output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=50)
pol,sfs,div = parseSfs(param=param,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=50)
sfs = convert.(Int64,cumulativeSfs(sfs))


Plots.theme(:vibrant)
p = hcat(α,α_nopos,a,b,alxSummStat)
lb=["Analytical α" "Analytical α_nopos" "Simulated α" "Simulated α_nopos" "Sampled α"]
l = @layout [a; b; c]
fn = plot(p,legend = :outertopright,label=lb)
fn1 = plot(α,label="Analytical α",legend = :outertopright)
fn2 = plot(alxSummStat,label="Sampled α",legend = :outertopright,linecolor="#d43625")
plot(fn,fn1,fn2,layout=l)

# x,y = analyticalAlpha(param=param)

function analyticalAlpha(;param::parameters)

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
    dnW = fNeg + fPosL
    dnS = fNeg + fPosH

	## Polymorphism
	neut = DiscSFSNeutDown(param)

	selH = DiscSFSSelPosDown(param,param.gH,param.pposH)
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL)
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL)
	splitColumns(matrix) = (view(matrix, :, i) for i in 1:size(matrix, 2))
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN))

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN
	selW = selL + selN
	selS = selH + selN

	## Outputs
	α = @. 1 - ((ds/dn) * (sel/neut))

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

	## Polymorphism
    sel_nopos = selN

    αW = param.alLow/param.alTot
	α_nopos  = @. 1 - (ds_nopos/dn_nopos) * (sel_nopos/neut)
	αW_nopos = α_nopos * αW
	αS_nopos  =  α_nopos - αW_nopos

	##########
	# Output #
	##########

	return (α,α_nopos)
	# return (α,α_nopos,[α[end] asymp1[1] c1[1] c2[1] c3[1]],[α_nopos[end] asymp2[1] c1[2] c2[2] c3[2]])
end


function alphaByFrequencies(param::parameters,divergence::Array{Int64,1},sfs::Array{Int64,2},bins::Int64,cutoff::Float64)

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
	neut = DiscSFSNeutDown(param)

	selH = DiscSFSSelPosDown(param,param.gH,param.pposH);
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL);
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL);

	sel = (selH+selL)+selN;

	cumulativePs = cumulativeSfs(neut)[:,1]
	cumulativePn = cumulativeSfs(sel)[:,1]

	## Outputs
	α = sampledAlpha(d=divergence,afs=sfs,λdiv=hcat(ds,dn),λpol=hcat(cumulativePs,cumulativePn),expV=false,bins=bins)
	α = view(α,1:trunc(Int64,param.nn*cutoff),:)

	αS, expectedDn, expectedDs, expectedPn, expectedPs, summStat = sampledAlpha(d=divergence,afs=sfs,λdiv=hcat(ds,dn),λpol=hcat(cumulativePs,cumulativePn),expV=true,bins=bins)
	# d=divergence;afs=sfs;λdiv=hcat(ds,dn);λpol=hcat(cumulativePs,cumulativePn);expV=true;bins=bins
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
	cumulativePn_nopos = cumulativeSfs(sel_nopos)[:,1]

	## Outputs
	α_nopos = sampledAlpha(d=divergence,afs=sfs,λdiv=hcat(ds,dn),λpol=hcat(cumulativePs,cumulativePn_nopos),expV=false,bins=bins)
	# α_nopos = view(α_nopos,1:trunc(Int64,param.nn * cutoff),:)

	##########
	# Output #
	##########
	Dn,Ds,Pn,Ps = expectedDn,expectedDs,sum(view(expectedPn,1,:),dims=2),sum(view(expectedPs,1,:),dims=2)

	alphas = round.(hcat(α[end], α_nopos[end] .- α[end], α_nopos[end]),digits=5)
	alphas = repeat(alphas,outer=[2,1])

	expectedValues = hcat(DataFrame(alphas),DataFrame(hcat(Dn,Ds,Pn,Ps)),DataFrame(permutedims(summStat)),makeunique=true)

	return (α,α_nopos,expectedValues)
end



function setPpos()
 	sc          = pyimport("scipy.optimize")
	pposL,pposH = sc.fsolve(solvEqns,(0.0,0.0))

	if pposL < 0.0
	 	pposL = 0.0
	end
	if pposH < 0.0
		pposH = 0.0
   end
	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results

	param.pposL,param.pposH = pposL, pposH
end


function resampleByIndex(;param::parameters,bArr::BitArray{1},alpha::Array{Float64,1},afs::Array{Int64,2},pn::Array{Float64,1},ps::Array{Float64,1},Pn::Array{Int64,2},Ps::Array{Int64,2},Dn::Array{Int64,1},Ds::Array{Int64,1},stats::Array{Float64,2},bins::Int64,cutoff::Float64)

	while sum(bArr) < size(alpha,1)
		id = findall(x ->x == false, bArr)

		if size(id,1) < 2
			id = id[1]
		end

		tmpPn, tmpPs    = poissonPolymorphism(observedValues=afs[:,id],λps=ps,λpn=pn)

		Pn[:,id] = tmpPn; Ps[:,id] = tmpPs;
		reduceExpPn = view(reduceSfs(tmpPn,bins)',1:bins,:);
		reduceExpPs = view(reduceSfs(tmpPs,bins)',1:bins,:);

		stats[:,id] = @. 1 - ((Ds[id]/Dn[id])' * (reduceExpPn./reduceExpPs));
		stats = round.(stats,digits=5);

		tmpAlpha = @. 1 - ((Ds[id]/Dn[id])' * (tmpPn/tmpPs))
		alpha[id,:] =  @view tmpAlpha[trunc(Int64,param.nn*cutoff),:]
		bArr[id,:] = alpha[id,:] .> 0

	end

	return alpha,stats
end


function scipyFit(alphaValues::Array{Float64,1})

	x = collect(1:size(alphaValues,1))
	py"""
	import numpy as np
	# from numba import *
	from scipy import optimize
	# @njit(cache=True)
	def exp_model(x,a,b,c):
		return a + b*np.exp(-x*c)
	# @jit(cache=True)
	def test(x1,al):
		res = {}
		model = optimize.curve_fit(exp_model,x1, al, method='dogbox')

		res['a'] = model[0][0]
		res['b'] = model[0][1]
		res['c'] = model[0][2]

		# alpha for predicted model
		res['alpha'] = exp_model(x1[-1], res['a'], res['b'], res['c'])
		return(res['alpha'])
	"""

	# plot(x,alphaTrim)
	# plot!(x,asympModel(x,fitted.param),legend=:bottomleft)

	return py"test"(PyObject(x),PyObject(alphaValues[:,1]))

end

function alphaByFrequenciesSampled(param::parameters,divergence::Array,sfs::Array{Int64,2},bins::Int64,cutoff::Float64)

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
	neut = DiscSFSNeutDown(param)

	selH = DiscSFSSelPosDown(param,param.gH,param.pposH);
	selL = DiscSFSSelPosDown(param,param.gL,param.pposL);
	selN = DiscSFSSelNegDown(param,param.pposH+param.pposL);

	sel = (selH+selL)+selN;

	cumulativePs = cumulativeSfs(neut)[:,1]
	cumulativePn = cumulativeSfs(sel)[:,1]

	## Outputs
	# α = sampledAlpha(d=divergence,afs=sfs,λdiv=hcat(ds,dn),λpol=hcat(neut,sel),expV=false,bins=bins)
	# α = view(α,1:trunc(Int64,param.nn*cutoff),:)

	# α, expectedDn, expectedDs, expectedPn, expectedPs, summStat = sampledAlpha(d=divergence,afs=sfs,λdiv=hcat(ds,dn),λpol=hcat(cumulativePs,cumulativePn),expV=true,bins=bins)
	expDn, expDs    = poissonFixation(observedValues=divergence,λds=ds,λdn=dn);
	expPn, expPs    = poissonPolymorphism(observedValues=sfs,λps=cumulativePs,λpn=cumulativePn);
	reduceExpPn = view(permutedims(reduceSfs(expPn,bins)),1:bins,:);
	reduceExpPs = view(permutedims(reduceSfs(expPs,bins)),1:bins,:);

	## Alpha from expected values. Used as summary statistics
	summStat = @. 1 - ((expDs/expDn)' * (reduceExpPn./reduceExpPs));
	summStat = round.(summStat,digits=5)

	tmp = @. 1 - ((expDs/expDn)' * (expPn/expPs));
	α = tmp[trunc(Int64,param.nn*cutoff),:]

	boolArr = (α .> 0) .& (α .< param.alTot)

	α, summStat = resampleByIndex(param=param,bArr=boolArr,alpha=α,afs=sfs,pn=cumulativePn,ps=cumulativePs,Pn=expPn,Ps=expPs,Dn=expDn,Ds=expDs,stats=summStat,bins=bins,cutoff=cutoff)

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
	cumulativePn_nopos = cumulativeSfs(sel_nopos)[:,1]

	## Outputs
	expPn_nopos, expPs_nopos    = poissonPolymorphism(observedValues=sfs,λps=cumulativePs,λpn=cumulativePn_nopos)
	tmp_nopos = @. 1 - ((expDs/expDn)' * (expPn_nopos/expPs_nopos))
	α_nopos = tmp_nopos[trunc(Int64,param.nn*cutoff),:]

	alBoolArr = α_nopos .> α
	while sum(alBoolArr) < size(α_nopos,1)
		id = findall(x ->x == false, alBoolArr)

		if size(id,1) < 2
			id = id[1]
		end

		tmpPn, tmpPs = poissonPolymorphism(observedValues=sfs[:,id],λps=cumulativePs,λpn=cumulativePn_nopos)
		expPn_nopos[:,id] = tmpPn; expPs_nopos[:,id] = tmpPs;

		tmpAlpha = @. 1 - ((expDs[id]/expDn[id])' * (tmpPn/tmpPs))
		α_nopos[id,:] =  @view tmpAlpha[trunc(Int64,param.nn*cutoff),:]
		alBoolArr[id,:] = α_nopos[id,:] .> α[id,:]

	end

	##########
	# Output #
	##########
	Dn,Ds,Pn,Ps = expDn,expDs,sum(view(expPn,1,:),dims=2),sum(view(expPs,1,:),dims=2)

	alphas = round.(hcat(α, α_nopos .- α, α_nopos),digits=5)
	# alphas = repeat(alphas,outer=[2,1])

	expectedValues = hcat(DataFrame(alphas),DataFrame(hcat(Dn,Ds,Pn,Ps)),DataFrame(permutedims(summStat)),makeunique=true)

	return (α,α_nopos,expectedValues)
end

function posArblib(param::parameters,gammaValue::Int64,ppos::Float64,binom::SparseMatrixCSC{Float64,Int64})

	if ppos == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		redPlus = phiReduction(param,gammaValue)

		# Solving sfs
		NN2 = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = round.(xa1/(NN2),digits=10)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		gammaCorrected = gammaValue*param.B

		if isinf(exp(2*gammaCorrected))
			# Checked mpmath, BigFloat, DecFP.Dec128, Quadmath.Float128
			# gammaExp1 = exp(ArbFloat(gammaCorrected*2,bits=24))
			# gammaExp2 = exp(ArbFloat(gammaCorrected*-2,bits=24))
			
			gammaExp1 = exp(Arblib.Arb(gammaCorrected*2))
			gammaExp2 = exp(Arblib.Arb(gammaCorrected*-2))

		else
			gammaExp1 = exp(gammaCorrected*2)
			gammaExp2 = exp(gammaCorrected*-2)
		end
        
		function pS(i::Float64,g1::T,g2::T,ppos::Float64) where {T<:Union{Float64,Arb}}

			num1 = g1
			num2 = @. 1 - (g2 * exp(1-i))
			num = @. num1*num2
			den = @. (g1 -1) * i *(1-i)
		
			tmp = @. Float64(ppos*0.5*(num/den))
			return tmp
		end

		# Allocating outputs
		solvedPositiveSfs::Array{Float64,1} = (1.0/(NN2)) * (pS.(xa2,gammaExp1,gammaExp2,ppos))
		replace!(solvedPositiveSfs, NaN => 0.0)

		# subsetDict = get(param.bn,param.B,1)
		# out               = (param.thetaMidNeutral)*redPlus*0.745*(subsetDict*solvedPositiveSfs)
		out::Array{Float64,1} = (param.thetaMidNeutral)*redPlus*0.745*(binom*solvedPositiveSfs)
		# out = out[2:end-1]

	end

	return out
end

########################################################
################################
###    Summary statistics    ###
################################

function poissonFixation(;observedValues::Array, λds::Float64, λdn::Float64,replicas::Int64=1)

	ds = λds / (λds + λdn)
	dn = λdn / (λds + λdn)
	observedValues = repeat(observedValues,1,1,replicas)
	# poissonS  = (ds .* observedValues) .|> Poisson
	# poissonD  = (dn .* observedValues) .|> Poisson
	# sampledDs = rand.(poissonS,1)
	# sampledDn = rand.(poissonD,1)

	sampledDs  = PoissonRandom.pois_rand.(ds .* observedValues)
	sampledDn  = PoissonRandom.pois_rand.(dn .* observedValues)

	out = sampledDn, sampledDs
	return out
end

function poissonPolymorphism(;observedValues::Array, λps::Array{Float64,1}, λpn::Array{Float64,1},replicas::Int64=1)

	λ1 = similar(λps);λ2 = similar(λpn)

	observedValues = repeat(observedValues,1,1,replicas)
	sampledPs = similar(observedValues)
	sampledPn = similar(observedValues)

	# Neutral λ
	λ1 = @. λps / (λps + λpn)
	# Selected λ
	λ2 = @. λpn / (λps + λpn)

	sampledPs =  PoissonRandom.pois_rand.(observedValues .* λ1)
	sampledPn =  PoissonRandom.pois_rand.(observedValues .* λ2)

	return (sampledPn, sampledPs)	
end

function sampledAlpha(;param::parameters,d::Array,afs::Array,λdiv::Array{Float64,2},λpol::Array{Float64,2},replicas::Int64=1)

	pn = λpol[:,2]
	ps = λpol[:,1]

	## Outputs
	expDn, expDs    = poissonFixation(observedValues=d,λds=λdiv[1],λdn=λdiv[2],replicas=replicas)
	expPn, expPs    = poissonPolymorphism(observedValues=afs,λps=ps,λpn=pn,replicas=replicas)

	## Alpha from expected values. Used as summary statistics
	αS = @. round(1 - ((expDs/expDn) * (expPn/expPs)),digits=5)
	αS = reshape(αS,size(αS,1),size(αS,2)*size(αS,3))'

	return αS,expDn,expDs,expPn,expPs
end

function analyticalAlpha(;param::parameters)

	##############################################################
						# Solve the model  #
	##############################################################
	B = param.B

	setThetaF!(param)
	thetaF = param.thetaF
	# Solve the probabilities of fixations without background selection
	## First set non-bgs
	param.B = 0.999
	## Solve the mutation rate
	setThetaF!(param)
	## Solve the probabilities
	setPpos!(param)
	# Return to the original values
	param.thetaF = thetaF
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
	fN_nopos     = fN*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fNeg_nopos   = fNeg*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fPosL_nopos  = fPosL*(param.thetaMidNeutral/2.)*param.TE*param.NN
	fPosH_nopos  = fPosH*(param.thetaMidNeutral/2.)*param.TE*param.NN

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
end

function alphaByFrequencies(param::parameters,divergence::Array,sfs::Array,dac::Array{Int64,1},replicas::Int64=1)

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
	neut::Array{Float64,1} = DiscSFSNeutDown(param,param.bn[param.B])
	# neut = param.neut[param.B]

	selH::Array{Float64,1} = DiscSFSSelPosDown(param,param.gH,param.pposH,param.bn[param.B])
	selL::Array{Float64,1} = DiscSFSSelPosDown(param,param.gL,param.pposL,param.bn[param.B])
	selN::Array{Float64,1} = DiscSFSSelNegDown(param,param.pposH+param.pposL,param.bn[param.B])
	tmp = cumulativeSfs(hcat(neut,selH,selL,selN))
	splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));

	neut, selH, selL, selN = splitColumns(tmp)
	sel = (selH+selL)+selN
	
	## Outputs
	α = @. 1 - (ds/dn) * (sel/neut)
	# α = view(α,1:trunc(Int64,param.nn*cutoff),:)

	alxSummStat, expectedDn, expectedDs, expectedPn, expectedPs = sampledAlpha(param=param,d=divergence,afs=sfs[dac],λdiv=hcat(ds,dn),λpol=hcat(neut[dac],sel[dac]),replicas=replicas)

	##################################################################
	# Accounting for for neutral and deleterious alleles segregating #
	##################################################################
	## Fixation
	fN_nopos       = fN*(param.thetaMidNeutral/2.)*param.TE*param.NN
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
	α_nopos    = @. 1 - (ds_nopos/dn_nopos) * (sel_nopos/neut)
	#=αW_nopos   = @. 1 - (ds_nopos/dn_nopos) * ((selN + selL)/neut)
	αS_nopos   = @. 1 - (ds_nopos/dn_nopos) * ((selN + selH)/neut)=#
	
	# αW_nopos   = α_nopos * αW
	# αS_nopos   = α_nopos * (1 - αW)

	# α_nopos    = @. 1 - (ds_nopos/dn_nopos) * (sel_nopos/neut)[dac]
	# amk,ci,model = asympFit(α_nopos)
	# αW_nopos   = amk * αW
	# αS_nopos   = amk * (1 - αW)

	##########
	# Output #
	##########
	#=alphas = round.(hcat(α_nopos[(param.nn-1)] * αW , α_nopos[(param.nn-1)] * (1 - αW), α_nopos[(param.nn-1)]), digits=5)
	expectedValues = hcat(alphas,permutedims(alxSummStat))=#
	alphas = round.(hcat(α_nopos[(param.nn-1)] * αW , α_nopos[(param.nn-1)] * (1 - αW), α_nopos[(param.nn-1)]), digits=5)
	expectedValues = hcat(repeat(alphas,replicas),alxSummStat)

	return (α,α_nopos,expectedValues)
end

function writeSummaryStatistics(fileName::String,summStat)

	names = collect('a':'z')
	for i in 1:size(summStat,1)
		write(fileName * "_" * names[i] * ".tsv", DataFrame(summStat[i:i,:]), delim='\t', append=true)
	end
end
