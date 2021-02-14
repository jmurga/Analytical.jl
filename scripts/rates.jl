using Distributed
addprocs()
@everywhere using Analytical

####
adap = Analytical.parameters(N=1000,n=500,gamNeg=-457,gL=10,gH=500,bRange=permutedims(push!(collect(0.1:0.025:0.975),0.999)),dac = [2,4,5,10,20,50,200,500,700])

convolutedSamples = Analytical.binomialDict()
Analytical.binomOp!(adap,convolutedSamples.bn);

@time df = Analytical.ratesToStats(param = adap,convolutedSamples=convolutedSamples,gH=[600,500,400,300],gL=collect(5:10),gamNeg=[-200,-300,-400,-500,-1000],iterations = 10^2,shape=adap.al,output="/home/jmurga/ratesNew.jld2");
run(`rm /home/jmurga/ratesNew.jld2`)

@time df = Analytical.ratesToStats(param = adap,convolutedSamples=convolutedSamples,gH=[600,500,400,300],gL=collect(5:10),gamNeg=[-200,-400,-500,-1000],iterations = 10^6,shape=adap.al,output="/home/jmurga/ratesNew.jld2");


