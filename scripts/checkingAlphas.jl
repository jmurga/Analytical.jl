using Analytical, CSV, DataFrames, RCall;R"""library(ggplot2);library(data.table)"""
# Open empirical data
# using CSV, DataFrames, Plots
# param = parameters(N=1000,n=50,B=0.999,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=0.4,alLow=0.2);param.nn=101 ;binomOp(param)

#=struct onlyone <: AbstractMatrix{Bool}
    v::Bool
end
function Base.iterate(o::onlyone, state=1)
      state == 1 ? o.v : !o.v, state+1
end
Base.size(o::onlyone) = (1,typemax(Int))
Base.length(o::onlyone) = typemax(Int)
Base.getindex(o::onlyone,i) = i == 1 ? o.v : !o.v
Base.getindex(o::onlyone,i,j) = j == 1 ? o.v : !o.v
=#
function readSimulations(;analysis::String,b::Float64,dac::Array{Int64,1},output::String="/home/jmurga/mkt/202004/results/simulations/alphasModel",path::String="/home/jmurga/mkt/202004/rawData/simulations/")
	#=	analysis="noDemog_0.4_0.3_0.2"
		path="/home/jmurga/mkt/202004/rawData/simulations/"
		output = "/home/jmurga/mkt/202004/results/simulations/alphasModel"
	=#

	f = path * "/" * split(analysis,"_")[1] * "/" *analysis
	outPlot = output .* "/" .* analysis .* ".svg"

	tmp = parse.(Float64,split(f,"_")[2:end])
	α   = tmp[1]
	αW  = tmp[2]
	bgs = tmp[3]

	#=	if bgs == 0.999
		bgs = (bgs - 0.199)
	else 
		bgs = bgs - 0.1
	end=#
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
	alphaCumu        = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)[dac]
	#=alphaReduced = hcat(alpha',alphaCumu')=#

	# Estimating sampled alpha_x
	br = convert(Array,append!(collect(0.1:0.05:0.95),0.999)')
	param = Analytical.parameters(N=5000,n=500,B=b,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=α,alLow=αW,Lf=2*10^5,bRange = br)
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

	# x1,y1,z1 = Analytical.alphaByFrequencies(param,d,sum(sfs[:,1:2],dims=2),500,0.999)
	tmp = []
	for i in 1:5
		x2,y2,z2 = Analytical.alphaByFrequencies(param,d,sum(rSfsCumu[:,1:2],dims=2),999,0.999)
		push!(tmp,z2)
	end
	out = reduce(vcat,tmp)[:,4:end][:,dac]'

	df = R"""
		d1 = as.data.table($out);
		d1$f = seq(1,length($dac))
		d1$estimation = 'Analytical estimation'
		d1 = melt(d1, id.vars=c('f','estimation'))

		d2 = as.data.table($alphaCumu);
		d2$f = seq(1,length($dac))
		d2$estimation = 'Input alpha(x) cumulative'
		d2 = melt(d2, id.vars=c('f','estimation'))
		
		df = rbind(d1,d2)
		df$analysis = $analysis
		p =  ggplot(df) + geom_line(aes(f,value,color=estimation)) + scale_color_manual(values=c('gray','#ab2710')) + theme_bw()
		ggsave(p,file=$outPlot)
		df
	"""
	return(rcopy(df))
	# CSV.write("/home/jmurga/mkt/202004/results/simulations/alphasModel/alphas" * split(f,"/")[end] *".tsv", df, delim='\t';header=true)=#
end

df1 = readSimulations(analysis="noDemog_0.4_0.1_0.2",b=0.1,dac=collect(1:999))
df2 = readSimulations(analysis="noDemog_0.4_0.3_0.2",b=0.15,dac=collect(1:999))
df3 = readSimulations(analysis="noDemog_0.4_0.1_0.999",b=0.999,dac=collect(1:999))
df4 = readSimulations(analysis="noDemog_0.4_0.3_0.999",b=0.75,dac=collect(1:999))

tmp = vcat(df1,df2,df3,df4)

R"""
df = $tmp
df$analysis = factor(df$analysis,levels=c('noDemog_0.4_0.1_0.2','noDemog_0.4_0.3_0.2','noDemog_0.4_0.1_0.999','noDemog_0.4_0.3_0.999'))

p =  ggplot(df) + geom_line(aes(f,value,color=estimation)) + scale_color_manual(values=c('gray','#ab2710')) + theme_bw() + facet_wrap(.~analysis,ncol=2) + ylim(c(-0.5,0.4));p

ggsave(p,file='/home/jmurga/mkt/202004/results/simulations/alphasSimulations/noDemogModeled.svg')
"""