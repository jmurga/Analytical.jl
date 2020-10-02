using Distributed
addprocs()
@everywhere using Analytical, DataFrames, CSV
# Set up model
adap = Analytical.parameters(N=500,n=500,gam_neg=-457, gL=10,gH=500,Lf=10^5,B=0.999,alTot=0.4,alLow=0.2)
Analytical.binomOp!(adap);

sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.2_0.999/sfs.tsv")))
rSfs = convert.(Int64,Analytical.reduceSfs(Analytical.cumulativeSfs(sfs[:,2:end]),100))'
sCumu = Analytical.cumulativeSfs(sfs[:,2:end])

sfsPos   = sCumu[:,1] + sCumu[:,2]
sfsNopos = sCumu[:,4] + sCumu[:,2]

divergence = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.2_0.999/div.tsv")))
d = [convert(Int64,sum(divergence[1:2]))]

alpha = @. 1 - (divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2])

inputAbc = DataFrame(alpha')


@time df = Analytical.summaryStats(adap,3,d,sfsPos);
@time df = Analytical.summaryStats(adap,583,d,sfsPos);

CSV.write("/home/jmurga/mkt/202004/rawData/summStat/noDemog/noDemog_0.4_0.2_0.999/noDemog_0.4_0.2_0.999.tsv", df, delim='\t',header=false);



## Open empirical data
# path= "/home/jmurga/mkt/202004/rawData/";suffix="txt";
# files = path .* filter(x -> occursin(suffix,x), readdir(path))[1]

# pol,sfs,div = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=100)
# sfs = sfs[:,1]
# div = [div[1]]

# CSV.write("/home/jmurga/mkt/202004/rawData/summStat/noDemog/noDemog_0.4_0.2_0.999/sfsnoDemog.tsv", inputAbc, delim='\t',header=false);
#
# function summStats(param::Analytical.parameters,iterations::Int64,div::Array,sfs::Array)
#
#     fac       = rand(-2:0.05:2,iterations)
#     afac      = @. 0.184*(2^fac)
#     bfac      = @. 0.000402*(2^fac)
#
#     alTot     = rand(collect(0.01:0.01:0.6),iterations)
#     lfac      = rand(collect(0.1:0.1:0.9),iterations)
#     alLow     = @. round(alTot * lfac,digits=5)
#     nParam = [param for i in 1:iterations]
#     nDiv = [div for i in 1:iterations]
#     nSfs = [sfs for i in 1:iterations]
#
#     wp = CachingPool(workers())
#     b = pmap(bgsIter,wp,nParam,afac,bfac,alTot,alLow,nDiv,nSfs);
# 	return(b)
# 	#=return(reduce(vcat,b))=#
# end
#
# @everywhere function bgsIter(param::Analytical.parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,div::Array,sfs::Array)
#
#     r = Array{Float64}(undef, 17, 100 + 3	)
#     # r = zeros(1,103)
#     param.al = afac; param.be = bfac;
#     param.alLow = alLow; param.alTot = alTot;
#     iter = 1
# 	for j in param.bRange
#         param.B = j
#
#         Analytical.set_theta_f!(param)
#         theta_f = param.theta_f
#         param.B = 0.999
#         Analytical.set_theta_f!(param)
#         Analytical.setPpos!(param)
#         param.theta_f = theta_f
#         param.B = j
#         x,y,z = Analytical.alphaByFrequencies(param,div,sfs,100,0.9)
#
#         r[iter,:] = z
#         iter = iter + 1;
#         # println(z)
#     end
#     # return(reduce(vcat,r))
#     return(r)
# end

# function test(p,d,s)
#     r=[]
#     for i in 1:1000
#         push!(r,bgsIter(p,d,s))
#     end
#     return(r)
# end

# reduce(vcat,reduce(vcat,a))


# Analytical.summaryStatistics(output, z)
