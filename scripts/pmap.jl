
using Distributed
addprocs()
@everywhere using Analytical
# Set up model
adap = Analytical.parameters(N=1000,n=661,gam_neg=-457, gL=10,gH=500,Lf=2*10^5,B=0.999,alTot=0.4,alLow=0.2)
Analytical.binomOp(adap);
# analyticalApproach(adap)[1:1000,:]


## Open empirical data
path= "/home/jmurga/mkt/202004/rawData/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

pol,sfs,div = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/data",sfsColumns=[3,5],divColumns=[6,7],bins=100)
sfs = sfs[:,1]
div = [div[1]]


function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array)

    fac       = rand(-2:0.05:2,iter)
    afac      = @. 0.184*(2^fac)
    bfac      = @. 0.000402*(2^fac)
    
    alTot     = rand(collect(0.05:0.01:0.4),iter)
    lfac      = rand(collect(0.1:0.1:1),iter)
    alLow     = @. round(alTot * lfac,digits=2)
    nParam = [param for i in 1:iter]
    nDiv = [div for i in 1:iter]
    nSfs = [sfs for i in 1:iter]

    wp = CachingPool(workers())
    b = pmap(bgsIter,wp,nParam,afac,bfac,alTot,alLow,nDiv,nSfs)
	return(b)
	# return(reduce(vcat,b))
end

@everywhere function bgsIter(param::Analytical.parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,div::Array,sfs::Array)

    r = DataFrame[]
    param.al = afac; param.be = bfac; 
    param.alLow = alLow; param.alTot = alTot; 

	for j in param.bRange
        param.B = j

        Analytical.set_theta_f(param)
        theta_f = param.theta_f
        param.B = 0.999
        Analytical.set_theta_f(param)
        Analytical.setPpos(param)
        param.theta_f = theta_f
        param.B = j
        x,y,z = Analytical.alphaByFrequencies(param,div,sfs,100,0.9)
        push!(r,z)
        # println(z)
    end
    # return(reduce(vcat,r))
    return(r)
end

@time df = summStats(adap,10,div,sfs);



# function test(p,d,s)
#     r=[]
#     for i in 1:1000
#         push!(r,bgsIter(p,d,s))
#     end
#     return(r)
# end

# reduce(vcat,reduce(vcat,a))


# Analytical.summaryStatistics(output, z)
