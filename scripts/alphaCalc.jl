using Analytical
using PyPlot

function asympPlot(alphaValues,originalAsymp)
	fig    = plt.figure(figsize=(8, 6))
	plt.style.use("seaborn")

	x      = 1:size(alphaValues)[1]
	y1     = alphaValues[:,1]
	y2     = alphaValues[:,2]

	plt.plot(x, y1,"-",color="#91d5db")
	plt.plot(x, y1,"o",color="#91d5db")
	plt.plot(x, y2,"-",color="#fab4be")
	plt.plot(x, y2,"o",color="#fab4be")
	plt.plot(x, fill(originalAsymp,size(alphaValues)[1]),"-.",color="gray")
	plt.grid()
	plt.xscale("log")
	plt.xlim(0,60)
	# plt.ylim((y1.min() - (y1.min()*50/100)),0.25)
	plt.xticks([1,2,5,10,20,50],[1,2,5,10,20,50])
	plt.grid()

	return(fig)
end

###### Fig 1. ######
B  = 0.999
al = 0.2
weak = 0.2
# Change adap to desired values.
# Not declared options reset to default values
Analytical.changeParameters(gam_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=0.999,bRange=append!(collect(0.2:0.05:0.95),0.999),pposL=0.001,pposH=0,N=1000,n=661,Lf=10^6,rho=0.001,TE=5.0,convoluteBinomial=true)

empiricalValues = Analytical.parseSfs("/home/jmurga/mktest/data/mk_with_positions_BGS.txt",[3,5],[6,7])
# Adap contains all parameters needed to estimate alpha under N conditions
println(adap)



function test(iter)

	# println(i)
	for i in 1:iter
		
		
		for j in adap.bRange
		# for j in [0.999]
		
			Analytical.changeParameters(gam_neg=-rand(80:200),gL=rand(10:20),gH=rand(100:500),alLow=rand(collect(0.1:0.1:0.4)),alTot=rand(collect(0.1:0.1:0.4)),theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=j,bRange=bRange=append!(collect(0.2:0.05:0.95),0.999),pposL=0.001,pposH=0,N=1000,n=661,Lf=10^6,rho=0.001,TE=5.0,convoluteBinomial=false)

			Analytical.set_theta_f()
			theta_f = adap.theta_f
			adap.B = 0.999
			Analytical.set_theta_f()
			Analytical.setPpos()
			adap.theta_f = theta_f
			adap.B = j

			x,y,z= Analytical.alphaByFrequencies(adap.gL,adap.gH,adap.pposL,adap.pposH,empiricalValues,"af")
			CSV.write("/home/jmurga/allProteins.csv", DataFrame(z), delim='\t', append=true)
		end

	end
end

alphaAll = ABCreg(data="/home/jmurga/data2.gz",prior="/home/jmurga/prior.gz", nparams=5, nsummaries=4, outputPath="/home/jmurga/", outputPrefix="newOut", tolerance=0.001, regressionMode="T",regPath="/home/jmurga/ABCreg/src/reg")

alpha = hcat(pos,nopos)
asympPlot(alpha,0.2)

###### Fig 2. ######
function bgsPlot(bgs1,bgs2,bgs3,bgs4,trueAlpha)
	fig    = plt.figure(figsize=(8, 6))
	plt.style.use("seaborn")

	x      = 1:length(bgs1)
	y1     = bgs1
	y2     = bgs2
	y3     = bgs3
	y4     = bgs4

	plt.plot(x, y1,"-",color="blue")
	plt.plot(x, y2,"-",color="orange")
	plt.plot(x, y3,"-",color="yellow")
	plt.plot(x, y4,"-",color="red")
	plt.plot(x, fill(trueAlpha,length(bgs1)),"-.",color="gray")
	plt.grid()
	plt.xscale("log")
	plt.xticks([1,2,5,10,20,50],[1,2,5,10,20,50])
	plt.grid()

	plt.xlim(0,60)

	return fig
end

###########################
B  = [0.4,0.6,0.8,0.999]
al = float(0.2)

###########################

function test2(alphaLow)
	results = Array{Array}(undef,length(B))
	for i in 1:length(B)

		println(B[i])
		Analytical.changeParameters(gam_neg=-83,gL=10,gH=500,alLow=alphaLow,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=B[i],pposL=0.001,pposH=0,N=500,n=25,Lf=10^6,L_mid=501,rho=0.001,al2= 0.0415,be2=0.00515625,TE=5.0,ABC=false)

		# here the software calculates the mutation rates corresponding to the desired terms
		Analytical.set_theta_f()
		theta_f = adap.theta_f
		adap.B = 0.999
		Analytical.set_theta_f()
		Analytical.setPpos()
		adap.theta_f = theta_f
		adap.B = B[i]
		# run the calucation

		x = Analytical.alphaByFrequencies(adap.gL,adap.gH,adap.pposL,adap.pposH,"nopos")
		results[i] = x
	end
	return results
end

strong = test2(0.0)
bgsPlot(strong[1],strong[2],strong[3],strong[4],0.2)
weak = test2(0.2)
bgsPlot(weak[1],weak[2],weak[3],weak[4],0.2)
