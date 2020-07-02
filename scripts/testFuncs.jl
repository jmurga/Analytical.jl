# Open empirical data
param = parameters(N=1000,n=661,B=0.999,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=0.8,alLow=0.4) ;binomOp(param)
# param = parameters(N=1000,n=50,B=0.999,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=0.4,alLow=0.2);param.nn=101 ;binomOp(param)
# path= "/home/jmurga/mktest/data/";suffix="txt";
# files = path .* filter(x -> occursin(suffix,x), readdir(path))
j = param.B
set_theta_f(param)
theta_f = param.theta_f
param.B = 0.999
set_theta_f(param)
setPpos(param)
param.theta_f = theta_f
param.B = j

# pol,sfs,div = parseSfs(param=param,data=files[1],output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=50)
pol,sfs,div = parseSfs(param=param,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=50)
sfs = convert.(Int64,cumulativeSfs(sfs))

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
