println("=====> Loading modules")
using Analytical
using CSV
import DataFrames.DataFrame

sampleSize = parse(Int64,ARGS[1])
estimations = parse(Int64,ARGS[2])
files = convert(String,ARGS[3]) .* readdir(ARGS[3])
output = ARGS[4]

println("=====> Loading model")

Analytical.changeParameters(N=1000,n=sampleSize,convoluteBinomial=true)
empiricalData = Analytical.parseSfs(files,[3,5],[6,7])

function solveEstimations(iter,data,out)
	for i in 1:iter		
		# for j in adap.bRange
		for j in [0.999]
		
			Analytical.changeParameters(gam_neg=-rand(80:200),gL=rand(10:20),gH=rand(100:500),alLow=rand(collect(0.1:0.1:0.4)),alTot=rand(collect(0.1:0.1:0.4)),theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=j,bRange=bRange=append!(collect(0.2:0.05:0.95),0.999),pposL=0.001,pposH=0,N=1000,n=sampleSize,Lf=10^6,rho=0.001,TE=5.0,convoluteBinomial=false)

			Analytical.set_theta_f()
			theta_f = adap.theta_f
			adap.B = 0.999
			Analytical.set_theta_f()
			Analytical.setPpos()
			adap.theta_f = theta_f
			adap.B = j

			x,y,z= Analytical.alphaByFrequencies(adap.gL,adap.gH,adap.pposL,adap.pposH,data,"af")
			CSV.write(out, DataFrame(z), delim='\t', append=true)
			# println(z)
		end
	end
end

println("=====> Solving alpha")

@time solveEstimations(estimations,empiricalData,output)