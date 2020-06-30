using Analytical, ProgressMeter
# Set up model
adap = Analytical.parameters(N=500,n=661,gam_neg=-457, gL=10,gH=500,Lf=2*10^5,B=0.999,alTot=0.4,alLow=0.4)
Analytical.binomOp(adap)
# analyticalApproach(adap)[1:1000,:]


## Open empirical data
path= "/home/jmurga/mkt/202004/rawData/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

pol,sfs,div = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=100)

sfs = convert.(Int64,Analytical.cumulativeSfs(sfs))

function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
	# @threads
	@showprogress for i in 1:iter
	# for i in 1:iter


		fac       = rand(-2:0.05:2)
		afac      = 0.184*(2^fac)
		bfac      = 0.000402*(2^fac)
		
		alTot     = rand(collect(0.05:0.01:0.4))
		alLow     = rand(collect(0.0:0.01:alTot))

		# println((thread=Threads.threadid(), iteration=i))
		bgsIter(param,afac,bfac,alTot,alLow,div,sfs,output,b,c)
	end
end

function bgsIter(param::Analytical.parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)

	for j in param.bRange
			# j = 0.999
			param.al = afac; param.be = bfac; 
			param.alLow = alLow; param.alTot = alTot; param.B = j

			Analytical.set_theta_f(param)
			theta_f = param.theta_f
			param.B = 0.999
			Analytical.set_theta_f(param)
			Analytical.setPpos(param)
			param.theta_f = theta_f
			param.B = j
			# x,y = Analytical.analyticalAlpha(param=param)
			x,y,z = Analytical.alphaByFrequencies(param,div,sfs,b,c)
			
			
			# println(z)

			Analytical.summaryStatistics(output, z)

		end
end
	

summStats(adap,1,div,sfs,"/home/jmurga/test",100,0.9)
summStats(adap,7400,div,sfs,"/home/jmurga/test",100,0.9)
# Custom function to perform 10^6 random solutions
function analyticalApproach(param)

	Analytical.set_theta_f(param)
	theta_f = param.theta_f
	println(theta_f)
	param.B = 0.999
	Analytical.set_theta_f(param)
	Analytical.setPpos(param)
	param.theta_f = theta_f
	adap.B = param.B
	x,y = Analytical.analyticalAlpha(param=param)

	return hcat(x,y)
end

a,b=analyticalApproach(100)

dfA = DataFrame(a[a[:,2] .!= 0,:],[:lastValue, :asymp, :c1,:c2,:c3]) |> stack
dfB = DataFrame(b[xb[:,2] .!= 0,:],[:lastValue, :asymp, :c1,:c2,:c3]) |> stack

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

# plotPosterior(posteriusing Analytical, ProgressMeter

# Set up model
adap = Analytical.parameters(N=500,n=661,gam_neg=-457, gL=10,gH=500)
Analytical.binomOp(adap)

# # Open empirical data
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

div,sfs,pol = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=50)

function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
	# @threads
	@showprogress for i in 1:iter
	# for i in 1:iter


		fac       = rand(-2:0.05:2)
		afac      = 0.184*(2^fac)
		bfac      = 0.000402*(2^fac)

		
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
			Analytical.setPpos(param)
			param.theta_f = theta_f
			param.B = j
			# x,y = Analytical.analyticalAlpha(param=param)
			x,y,z = Analytical.alphaByFrequencies(param,div,sfs,b,c)
			
			
			# println(z)
			if sum(convert(Array,z[1:1,1:3]) .> 0) < 3
				# println("param.B=" * string(param.B) * ";param.al=" *string(param.al) *";param.be=" *string(param.be) * ";param.alTot="* string(param.alTot) *";param.alLow=" *string(param.alLow))
				# println(hcat(x[end],y[end]))
				# println(param)
				continue
			else
				Analytical.summaryStatistics(output, z)
			end

		end
	end
end
orVip[1],(600,400))
