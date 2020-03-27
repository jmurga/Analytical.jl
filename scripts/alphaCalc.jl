using Analytical
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


B  = 0.999
al = 0.2
weak = 0.2

# Adap contains all parameters needed to estimate alpha under N conditions
println(adap)

# Change adap to desired values.
# Not declared options reset to default values
changeParameters(gam_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=0.999,pposL=0.001,pposH=0,N=500,n=25,Lf=10^6,L_mid=501,rho=0.001,al2= 0.0415,be2=0.00515625,TE=5.0,ABC=false)


set_theta_f()
theta_f = adap.theta_f
adap.B = 0.999
set_theta_f()
setPpos()
adap.theta_f = theta_f
adap.B = B
# run the calucation
pos,nopos = alphaByFrequencies(adap.gL,adap.gH,adap.pposL,adap.pposH)

alpha = hcat(pos,nopos)
