using QuadGK
using SpecialFunctions
using Roots

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


function eMKT(sfs, div, m,cutoff):

	p_0 = convert(Integer,sfs[:,3] |> sum)
    p_i = convert(Integer,sfs[:,2] |> sum)
    d_0 = convert(Integer,div[2])
    d_i = convert(Integer,div[1])
    m_0 = convert(Integer,m[2])
 	m_i = convert(Integer,m[1])

    # divergence metrics
    ka = di / mi
    ks = d0 / m0
    omega = ka/ks

    ### Estimating alpha with pi/p0 ratio
    piMinus   = sfs[sfs[:,1] .<= cutoff,2] |> sum
    piGreater = sfs[sfs[:,1] .> cutoff,2] |> sum
    p0Minus   = sfs[sfs[:,1] .<= cutoff,3] |> sum
    p0Greater = sfs[sfs[:,1] .> cutoff,3] |> sum

	ratiop0 = p0Minus / p0Greater
	
    deleterious = piMinus - (piGreater * ratiop0)
    piNeutral   = convert(Integer,round(p_i - deleterious))

    alpha = 1 - (((p_i - deleterious) / p_0) * (d_0 / d_i))

    ## Estimation of b: weakly deleterious
    b = (deleterious / p_0) * (m0 / mi)

    ## Estimation of f: neutral sites
    f = (m0 * piNeutral) / (mi * p_0)

    ## Estimation of d, strongly deleterious sites
    d = 1 - (f + b)

    # pvalue =  HypothesisTests.FisherExactTest(p_0, d_0, piNeutral, d_i)


    return alpha

