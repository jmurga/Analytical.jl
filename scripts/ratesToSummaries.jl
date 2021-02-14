using Distributed
addprocs()
@everywhere using Analytical

h5file = jldopen(ratesFile)

adap = Analytical.parameters(N=1000,n=500)
adap.dac = h5file["accesToDac"]

sfs        = Array(CSV.read(sFiles,DataFrame))
divergence = Array(CSV.read(dFiles,DataFrame))

scumu = Analytical.cumulativeSfs(sfs)
s = sum(scumu[:,2:3],dims=2)[adap.dac]
d = [sum(divergence[1:2])]


summstat = Analytical.summaryStatsFromRates(param = adap,rates=h5file,divergence=[d],sfs=[s],summstatSize=10^5,replicas=1)