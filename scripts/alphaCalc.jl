using Analytical

# Set up sfs
Analytical.changeParameters(N=1000,n=661,diploid=true,convoluteBinomial=true)

# Open empirical data
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

empiricalValues = Analytical.parseSfs(data=files,output="data.tsv",sfsColumns=[3,5],divColumns=[6,7])

# Custom function to perform 10^6 random solutions
function summStats(iter,data,output)
	# @threads
	for i in 1:iter
		
		gam_neg=-rand(80:400)
		gL=rand(10:20)
		gH=rand(200:500)
		alLow=rand(collect(0.0:0.1:0.2))
		alTot=rand(collect(0.0:0.2:0.2))

		@threads for j in adap.bRange
			Analytical.changeParameters(gam_neg=gam_neg,gL=gL,gH=gH,alLow=alLow,alTot=alTot,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=j,bRange=adap.bRange,pposL=0.001,pposH=0.0,N=1000,n=661,Lf=10^6,rho=adap.rho,TE=5.0,convoluteBinomial=false)

			Analytical.set_theta_f()
			theta_f = adap.theta_f
			adap.B = 0.999
			Analytical.set_theta_f()
			Analytical.setPpos()
			adap.theta_f = theta_f
			adap.B = j

			x,y,z= Analytical.alphaByFrequencies(gammaL=adap.gL,gammaH=adap.gH,pposL=adap.pposL,pposH=adap.pposH,observedData=data)
			Analytical.summaryStatistics(output, z)
			
		end
	end
end

@time summStats(58824,empiricalValues,"/home/jmurga/prior.tsv")
