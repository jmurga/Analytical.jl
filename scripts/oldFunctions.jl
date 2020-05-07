
################################
######## Old functions  ########
################################
# function changeParameters(;γ_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,bRange=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.999],B=0.95,pposL=0.001,pposH=0,N=500,n=25,Lf=10^6,L_mid=501,ρ=0.001,al2= 0.0415,be2=0.00515625,TE=5.0,ABC=false,convoluteBinomial=true)
#
# 	adap.gam_neg           = γ_neg
# 	adap.gL                = gL
# 	adap.gH                = gH
# 	adap.alLow             = alLow
# 	adap.alTot             = alTot
# 	adap.theta_f           = theta_f
# 	adap.theta_mid_neutral = theta_mid_neutral
# 	adap.al                = al
# 	adap.be                = be
# 	adap.bRange            = bRange
# 	adap.B                 = B
# 	adap.pposL             = pposL
# 	adap.pposH             = pposH
# 	adap.N                 = N
# 	adap.n                 = n
# 	adap.Lf                = Lf
# 	adap.L_mid             = L_mid
# 	adap.rho               = ρ
# 	adap.al2               = al2
# 	adap.be2               = be2
# 	adap.TE                = TE
# 	adap.ABC               = ABC
#
# 	adap.NN = N*2
# 	adap.nn = n*2
#
# 	if convoluteBinomial == true
# 		adap.bn = Dict(bRange[i] => binomOp(bRange[i]) for i in 1:length(bRange))
# 	end
# end
# function setPpos_nlsolve()
#
# 	function f!(F,x)
# 		F[1] = alphaExpSimTot(x[1],x[2])-adap.alTot
# 		F[2] = alphaExpSimLow(x[1],x[2])-adap.alLow
# 	end
#
# 	# Scipy probably cannot solve due to floats, Julia does so I implemented the same version forcing from the original results
# 	pposL,pposH = NLsolve.nlsolve(f!,[ 0.00001; 0.000001]).zero
#
# 	if pposL < 0.0
# 	 	pposL = 0.0
# 	end
# 	if (pposH < 0.0 || pposH < 9e-15)
# 		pposH = 0.0
# 	end
# 	adap.pposL,adap.pposH = pposL, pposH
# end

# function set_Lf()
#
# 	i(L)   = Br(L,theta_f,adap,adap)-B
# 	tmpLf  = Roots.find_zero(i,100)
# 	Lf     = convert(Int64,round(tmpLf))
#
# 	return Lf
# end

# function summaryStatistics(fileName,simulationName,alphaPos,alphaNopos)
#
# 	if isfile(fileName)
# 		jldopen(fileName, "a+") do file
# 			group = JLD2.Group(file, simulationName*string(rand(Int64)))
# 			group["alphaPos"] = alphaPos
# 			group["alphaNopos"] = alphaNopos
# 		end
# 	else
# 		jldopen(fileName, "w") do file
# 			group = JLD2.Group(file, simulationName*string(rand(Int64)))
# 			group["alphaPos"] = alphaPos
# 			group["alphaNopos"] = alphaNopos
# 		end
# 	end
# end

# function alphaByFrequenciesNoPositive(gammaL,gammaH,pposL,pposH)
#
# 	fN = adap.B*fixNeut()*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
# 	fNeg = adap.B*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN*fixNegB(0.5*pposH+0.5*pposL)
# 	fPosL = fixPosSim(gammaL,0.5*pposL)*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
# 	fPosH = fixPosSim(gammaH,0.5*pposH)*(adap.theta_mid_neutral/2.)*adap.TE*adap.NN
#
# 	neut = cumulativeSfs(DiscSFSNeutDown())
# 	selN = cumulativeSfs(DiscSFSSelNegDown(pposL+pposH))
#
# 	# Outputs
# 	ret = Array{Float64}(undef, adap.nn - 1)
# 	sel = Array{Float64}(undef, adap.nn - 1)
# 	for i in 1:length(ret)
# 		sel[i] = (selH[i]+selL[i])+selN[i]
# 		ret[i] = float(1.0 - (fNNopos/(fPosLNopos + fPosHNopos+  fNegNopos+0.0))* selNopos[i]/neutNopos[i])
# 	end
#
# 	return ret
# end

# set_theta_f()
# setPpos()

# adap.bn = Analytical.binomOp()
# function test()
# 		Analytical.set_theta_f()
# 		theta_f = adap.theta_f
# 		adap.B = 0.999
# 		Analytical.set_theta_f()
# 		Analytical.setPpos()
# 		adap.theta_f = theta_f
# 		adap.B = 0.4
# 		c,d,e,f,g,h,i,j = Analytical.alphaByFrequencies(adap.gL,adap.gH,adap.pposH,adap.pposL)
# end


#
# alphaExpSimLow(pposL,pposH) = fixPosSim(adap.gL,0.5*pposL)/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH) + fixNegB(0.5*pposL+0.5*pposH))
#
# alphaExpSimTot(pposL,pposH) = (fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH))/(fixPosSim(adap.gL,0.5*pposL)+fixPosSim(adap.gH,0.5*pposH)+fixNegB(0.5*pposL+0.5*pposH))

# function binomOp()
#
# 	sc = pyimport("scipy.stats")
# 	NN2 = convert(Int64, round(adap.NN*adap.B, digits=0))
#
# 	samples =  [i for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = [j for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = samplesFreqs./NN2
#
# 	return sc.binom.pmf(samples,adap.nn,samplesFreqs)
# 	NN2          = convert(Int64, round(adap.NN*adap.B, digits=0))
# 	samples      =  [i for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = [j for i in 0:adap.nn, j in 0:NN2]
# 	samplesFreqs = samplesFreqs/NN2
#
# 	f(x) = Distributions.Binomial(adap.nn,x)
# 	z    = samplesFreqs .|> f
#
# 	out  = Array{Float64}(undef,(adap.nn+1,adap.NN))
# 	out  = Distributions.pdf.(z,samples)
#
# 	return out
# end
