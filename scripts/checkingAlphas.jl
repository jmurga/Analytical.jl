# Open empirical data
using Analytical, CSV, DataFrames, Plots
# param = parameters(N=1000,n=50,B=0.999,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=0.4,alLow=0.2);param.nn=101 ;binomOp(param)

struct onlyone <: AbstractMatrix{Bool}
    v::Bool
end
function Base.iterate(o::onlyone, state=1)
      state == 1 ? o.v : !o.v, state+1
end
Base.size(o::onlyone) = (1,typemax(Int))
Base.length(o::onlyone) = typemax(Int)
Base.getindex(o::onlyone,i) = i == 1 ? o.v : !o.v
Base.getindex(o::onlyone,i,j) = j == 1 ? o.v : !o.v

function readSimulations(model::String,path::String="/home/jmurga/mkt/202004/rawData/simulations/")

	dir = path * model

	files = dir .* "/" .* readdir(dir)
	splitValues = readdir(dir)

	for f in files
		println(split(f,"/")[end])
		tmp = parse.(Float64,split(f,"_")[2:end])
		α   = tmp[1]
		αW  = tmp[2]
		bgs = tmp[3]

		# Open files and make inputs
		sfs = convert(Array,DataFrame!(CSV.File(f * "/sfs.tsv")))
		divergence = convert(Array,DataFrame!(CSV.File(f * "/div.tsv")))

		sfs = sfs[:,2:end]
		sCumu = convert.(Int64,Analytical.cumulativeSfs(sfs))

		sfsPos   = sCumu[:,1] + sCumu[:,2]
		sfsNopos = sCumu[:,4] + sCumu[:,2]
		d = [convert(Int64,sum(divergence[1:2]))]

		# Estimating input alpha_x
		rSfs         = Analytical.reduceSfs(sfs,999)'
		rSfsCumu     = Analytical.cumulativeSfs(Analytical.reduceSfs(sfs,999)')
		alpha        = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'
		alphaCumu        = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)'
		alphaReduced = hcat(alpha',alphaCumu')

		# Estimating sampled alpha_x
		param = Analytical.parameters(N=500,n=500,B=0.3,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=α,alLow=αW,Lf=2*10^5)
		# param.nn = param.nn+1
		Analytical.binomOp!(param)

		j = param.B
		Analytical.set_theta_f!(param)
		theta_f = param.theta_f
		param.B = 0.999
		Analytical.set_theta_f!(param)
		Analytical.setPpos!(param)
		param.theta_f = theta_f
		param.B = j

		x1,y1,z1 = Analytical.alphaByFrequencies(param,d,sum(sfs[:,1:2],dims=2),999,0.999)
		x2,y2,z2 = Analytical.alphaByFrequencies(param,d,sum(sCumu[:,1:2],dims=2),999,0.999)
		asymp1 = Analytical.asympFit(alpha')[2]

		# Df to plot
		dfSampled = hcat(alphaReduced,z1[:,4:end]',z2[:,4:end]')
		p1        = plot(dfSampled,label = ["Input α(x)" "Input α(x) cumulative" "sampled α(x)" "sampled α(x) cumulative"],title = "noDemog: " * "α=" * string(α) * ";αW=" * string(αW) * ";B="* string(bgs),linewidth=1.5,legend=:outerright,ylim=(-1,0.5))
		p1        =hline!([asymp1],label="Asymptotic α",linecolor =:black,linestyle =:dot)

		out1 = []
		# out2 = []
		for i in 1:50
			# x1,y1,z1 = Analytical.alphaByFrequencies(param,d,sum(sfs[:,1:2],dims=2),100,0.999)
			x1,y1,z1 = Analytical.alphaByFrequencies(param,d,sum(sCumu[:,1:2],dims=2),999,0.999)
			# push!(out2,z2[4:end])
			push!(out1,z1[4:end])
		end
		samplingDf = hcat(alphaCumu',reduce(hcat,out1))
		p2 = plot(samplingDf[:,2:end];label="Sampled α(x) x 50", linecolor =:gray,linewidth=0.2, primary=onlyone(true))
		p2 = plot!(samplingDf[:,1];label="Input α(x)",legend=:outerright,linewidth=2,linecolor=:black,ylim=(-1,0.5))


		l = @layout [a; b]
		plot(p1,p2,layout=l)


		savefig("/home/jmurga/mkt/202004/results/simulations/alphasModel/" * split(f,"/")[end] *".svg")
		df = DataFrame(dfSampled)
		rename!(df,[:input_α,:input_α_cumu,:sampled_α,:sampled_α_cumu])

		CSV.write("/home/jmurga/mkt/202004/results/simulations/alphasModel/alphas" * split(f,"/")[end] *".tsv", df, delim='\t';header=true)
	end
end
