using Analytical, ProgressMeter

# Set up model
adap = Analytical.parameters(N=1000,n=661,gam_neg=-457, gL=10,gH=500)
Analytical.binomOp(adap)

# # Open empirical data
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

empiricalValues = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=50)

adap.n = 25
adap.nn = 50
Analytical.binomOp(adap);
# # Custom function to perform 10^6 random solutions
function summStats(param::Analytical.parameters,iter::Int64,data::Array,output::String,b::Int64,c::Float64)
	# @threads
	@showprogress for i in 1:iter


		fac       = rand(-2:0.05:2)
		afac      = 0.184*(2^fac)
		# afac      = 0.184
		bfac      = 0.000402*(2^fac)
		# bfac      = 0.000402

		alTot     = rand(collect(0.05:0.01:0.4))
		alLow     = rand(collect(0.0:0.01:alTot))

		
		for j in param.bRange
		# j = 0.999

			param.al = afac; param.be = bfac; 
			param.alLow = alLow; param.alTot = alTot; param.B = j

			
			Analytical.set_theta_f(param)
			theta_f = param.theta_f
			param.B = 0.999
			Analytical.set_theta_f(param)
			Analytical.setPpos(param=param)
			param.theta_f = theta_f
			param.B = j
			x,y,z = Analytical.alphaByFrequencies(param,data,b,c)
			# x,y,z = Analytical.analyticalAlpha(param=param)
			
			if sum(convert(Array,z[1:1,1:3]) .> 0) < 3
				continue
			else
				Analytical.summaryStatistics(output, z)
			end

		end
	end
end

summStats(adap,1,empiricalValues,"/home/jmurga/test",50,0.98)
summStats(adap,117648,empiricalValues,"/home/jmurga/priorSample",20,0.999)


# Custom function to perform 10^6 random solutions
function analyticalApproach(iter::Int64)
	# @threads
	result = Array{Float64}(undef,iter*17,5)
	resultNopos = Array{Float64}(undef,iter*17,5)

	it = 1
	@showprogress for i in 1:iter
	# for i in 1:iter

		gam_neg   = -457
		gL        = 10
		gH        = 500

		fac       = rand(-2:0.05:2)
		afac      = 0.184*(2^fac)
		# afac      = 0.184
		bfac      = 0.000402*(2^fac)
		# bfac      = 0.000402

		alTot     = rand(collect(0.05:0.01:0.4))
		alLow     = rand(collect(0.0:0.01:alTot))


		for j in adap.bRange
		# j = 0.999

			Analytical.changeParameters(gam_neg=gam_neg,gL=gL,gH=gH,alLow=alLow,alTot=alTot,theta_f=1e-3,theta_mid_neutral=1e-3,al=afac,be=bfac,B=j,bRange=adap.bRange,pposL=0.001,pposH=0.0,N=1000,n=661,Lf=10^6,rho=adap.rho,TE=5.0,diploid=true,convoluteBinomial=false)

			Analytical.set_theta_f()
			theta_f = adap.theta_f
			adap.B = 0.999
			Analytical.set_theta_f()
			Analytical.setPpos()
			adap.theta_f = theta_f
			adap.B = j
			x1,y1,a1,a2 = Analytical.analyticalAlpha(gammaL=adap.gL,gammaH=adap.gH,pposL=adap.pposL,pposH=adap.pposH)

			result[it,:] = a1
			resultNopos[it,:] = a2
			# x,y,z = alphaByFrequencies(gammaL=adap.gL,gammaH=adap.gH,pposL=adap.pposL,pposH=adap.pposH,observedData=data,bins=20)

			it = it + 1
			# summaryStatistics(output, z)

		end

	end
	return result,resultNopos
end

a,b=analyticalApproach(100)

dfA = DataFrame(a[a[:,2] .!= 0,:],[:lastValue, :asymp, :c1,:c2,:c3]) |> stack
dfB = DataFrame(b[b[:,2] .!= 0,:],[:lastValue, :asymp, :c1,:c2,:c3]) |> stack

p2 = boxplot(["lastValue" "asymp" "90%" "80%" "75%"],b)


# using Plots, Plots.PlotMeasures, StatsPlots, LaTeXStrings, ProgressMeter, BenchmarkTools, DataFrames

# function plotPosterior(data,imgSize)

#     Plots.gr()
#     Plots.theme(:wong2)

#     p1 = StatsPlots.density(data,
#                             legend = :topleft,
#                             fill=(0, 0.75),
#                             linecolor=["#30504f" "#e2bd9a" "#ab2710"],
#                             #linecolor=["#E1B16A" "#31A2AC" "#F62A00"],
#                             fillcolor=["#30504f" "#e2bd9a" "#ab2710"],
#                             xlabel = L"\alpha",
#                             label = [L"\alpha_S" L"\alpha_W" L"\alpha"],
#                             ylabel = "Posterior density",
#                             lw = 2,
#                             fmt = :svg,
#                             size=imgSize,
#                         )
#     return p1
# end

# posteriorAll, resultAll = Analytical.ABCreg(data="/home/jmurga/dataAbc/dataAll.gz",prior="/home/jmurga/dataAbc/priorAll.gz", nparams=3, nsummaries=24, outputPath="/home/jmurga/dataAbc/", outputPrefix="all", tolerance=0.00025, regressionMode="T",regPath="/home/jmurga/ABCreg/src/reg")

# plotPosterior(posteriorAll[1],(600,400))

# posteriorVip, resultVip = Analytical.ABCreg(data="/home/jmurga/dataAbc/dataVip.gz",prior="/home/jmurga/dataAbc/priorVip.gz", nparams=3, nsummaries=24, outputPath="/home/jmurga/dataAbc/", outputPrefix="vip", tolerance=0.001, regressionMode="T",regPath="/home/jmurga/ABCreg/src/reg")

# plotPosterior(posteriorVip[1],(600,400))
