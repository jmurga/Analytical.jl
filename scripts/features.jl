using QuadGK
using SpecialFunctions
using Roots
using RCall

function Gammadist(gamma)
	return ((adap.be^adap.al)/SpecialFunctions.gamma(adap.al))*(gamma^(adap.al-1))*exp(-adap.be*gamma)
end

function pip0(gamma)
	U = 4*adap.theta_f*adap.Lf/(2.0*adap.NN)
	R = 2*adap.Lf*adap.rho/(2.0*adap.NN)
	return Gammadist(gamma)*exp(-(Gammadist(gamma)*U/(2.0*adap.NN))/(gamma/(adap.NN+0.0)+R/(2.0*adap.NN)))
end

function intpip0()

	f(gam) = pip0(gam)
	return QuadGK.quadgk(f,0,1000)[1]
end

function calcBGam(L,alpha,beta,theta)
	# not sure whats going on here, seems to overestimate or under sometimes
	u = 2*theta/(2.0*adap.NN)
	r = adap.rho/(2.0*adap.NN)

	a = -1.0*(2^(alpha+1))*exp(L*adap.NN*r*beta)*L*adap.N*((1.0/(L*adap.N*r))^(1-alpha))*u*(beta^alpha)
	b =  convert(Float64,SpecialFunctions.gamma_inc(1-alpha,L*adap.NN*r*beta,1)[2])
	#fudge = 1-gamInt.cdf(0.2,a=adap.al2,scale=1./adap.be2)
	#fudge = 1.
	fudge = 0.25

	c = exp(a*b*fudge)

	return c
end

function set_theta_f_gam()
	i(theta) = calcBGam(adap.Lf,adap.al2,adap.be2,theta)-adap.B
	# theta_f  = fsolve(lambda theta ,0.0000001)
	theta_f  = Roots.find_zero(i,0.0000001)
	adap.theta_f = theta_f
end

function calcB(L,theta)
	t = -1.0*adap.gam_neg/(adap.NN+0.)
	u = 2.0*theta/(2.0*adap.NN)
	r = adap.rho/(2.0*adap.NN)

	#return u*t/((t+r*L)^2)
	#return u*t/(2*(t+(1.-np.exp(-2*r*L))/2.)^2)

	# Nordborg models
	#####return u/(2*r*t*(1+((1.-np.exp(-2*r*L))/2.)*(1-t)/t)^2)
	#return u/(t*(1+((1-np.exp(-2.0*r*L))/2.)*(1-t)/t)^2)
	#return u/(t*(1+r*L*(1-t)/t)^2)
	#####return u/(t*(1+(np.exp(-2*r*L)/2)/t)^2)
end

# Br() needs L and theta, not working right now
function get_B_vals()
	ret = Array{Float64}(undef,30)
	for i in 20:50
		L = convert(Int64,round(1.3^i,digits=0))
		ret[i] = (Br(L),L)
	end
	return ret
end

"""
	Estimating and plotting MAP using locfit and ggplot2 in R. It assume your folder contains the posterior estimated through ABCreg
"""
function plotMap(;analysis::String,output::String)

	try
		@eval using RCall
		@eval R"""library(ggplot2);library(abc)"""

		out = filter(x -> occursin("post",x), readdir(analysis,join=true))
		out = filter(x -> !occursin(".1.",x),out)

		open(x) = Array(CSV.read(GZip.open(x),DataFrame,header=false))
		posteriors = open.(out)

		maxp = DataFrame(Array{Float64,2}(undef,size(posteriors,1),5),[:aw,:as,:a,:gamNeg,:shape])
		R"""getmap <- function(df){
				temp = as.data.frame(df)
			    d <-locfit(~temp[,1],temp);
			    map<-temp[,1][which.max(predict(d,newdata=temp))]
			}"""
		getmap(x) = rcopy(R"""matrix(apply($x,2,getmap),nrow=1)""")
		tmp = getmap.(posteriors)
		maxp = DataFrame(vcat(tmp...),[:aw,:as,:a,:gamNeg,:shape])
	
		al  = maxp[:,1:3]
		gam  = maxp[:,4:end]
		p = rcopy(R"""al = $al
			dal = melt(al)
			pal = ggplot(dal) + geom_density(aes(x=value,fill=variable),alpha=0.5) + scale_fill_manual(values=c('#30504f', '#e2bd9a', '#ab2710'))
			ggsave(pal,filename=paste0($output))
			""")
		return(maxp)
	catch
		println("Please install R, ggplot2 and abc in your system before execute this function")
	end
end