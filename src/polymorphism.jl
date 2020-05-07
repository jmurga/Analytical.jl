
################################
######    Polymorphism    ######
################################
# Expected number of polymorphism above frequency x from the standard diffusion theory
# f(x) = ∫(s) θs* (1/x(1-x)) * ( (ℯ^(Ns)*(1-ℯ^(-4Ns(1-x))) ) / (ℯ^(4Ns-1)))
# Convolote with the binomial to obtain the downsampled sfs
# E[P(x)] = ∑_x=x*→x=1 fB(x*)
############Neutral#############

"""
	Expected neutral allele frequency rate reduce by background selection. The spectrum depends on the number of individual: adap.NN.
	```math
		\\mathbb{E}\\left[P_{(x)}\\rigth] = \\sum{x^{*}=x}{x^{*}=1} f_{B}(x)
	```
# Return:
	- Array{Float64}(adap.NN · adap.B,1)
"""
function DiscSFSNeutDown()

	NN2 = convert(Int64,round(adap.NN*adap.B))

	# function neutralSfs(i)
	# 	if i > 0 && i < NN2
	# 		 return 1.0/(i)
	# 	end
	# 	return 0.0
	# end

	neutralSfs(i)    = 1.0/i

	x                = [convert(Float64,i) for i in 0:NN2]
	solvedNeutralSfs = x .|> neutralSfs
	out              = adap.B*(adap.theta_mid_neutral)*0.255*(adap.bn[adap.B]*solvedNeutralSfs)

	return 	out[2:lastindex(out)-1]
end

############Positive############
# Variable gamma in function changed to gammaValue to avoid problem with exported SpecialFunctions.gamma

function DiscSFSSelPosDown(gammaValue::Int64,ppos::Float64)

	if ppos == 0.0
		out = zeros(Float64,adap.nn + 1)
	else
		S        = abs(adap.gam_neg/(1.0*adap.NN))
		r        = adap.rho/(2.0*adap.NN)
		μ        = adap.theta_f/(2.0*adap.NN)
		s        = gammaValue/(adap.NN*1.0)
		Ψ0       = SpecialFunctions.polygamma(1,(s+S)/r)
		Ψ1       = SpecialFunctions.polygamma(1,1.0+(r*adap.Lf+s+S)/r)
		red_plus = ℯ^(-2.0*S*μ*(Ψ0-Ψ1)/(r^2))
				   
		# Solving sfs
		NN2 = convert(Int64,round(adap.NN*adap.B,digits=0))
		xa  = [i for i in 0:NN2]
		xa  = xa/(NN2+0.0)

		function positiveSfs(i,gammaCorrected=gammaValue*adap.B,ppos=ppos)
			if i > 0 && i < 1.0
				return ppos*0.5*(ℯ^(2*gammaCorrected)*(1-ℯ^(-2.0*gammaCorrected*(1.0-i)))/((ℯ^(2*gammaCorrected)-1.0)*i*(1.0-i)))
			end
			return 0.0
		end

		solvedPositiveSfs = (1.0/(NN2+0.0)) * (xa .|> positiveSfs)
		out               = (adap.theta_mid_neutral)*red_plus*0.745*(adap.bn[adap.B]*solvedPositiveSfs)
	end

	return out[2:lastindex(out)-1]
end

######Slightly deleterious######
function DiscSFSSelNegDown(ppos::Float64)
	out = adap.B*(adap.theta_mid_neutral)*0.745*(adap.bn[adap.B]*DiscSFSSelNeg(ppos))
	return out[2:lastindex(out)-1]
end

function DiscSFSSelNeg(ppos::Float64)

	beta     = adap.be/(1.0*adap.B)
	NN2      = convert(Int64, round(adap.NN*adap.B, digits=0))
	xa       = [round(i/(NN2+0.0),digits=6) for i in 0:NN2]

	z(x,ppos=ppos) = (1.0-ppos)*(2.0^-adap.al)*(beta^adap.al)*(-SpecialFunctions.zeta(adap.al,x+beta/2.0) + SpecialFunctions.zeta(adap.al,(2+beta)/2.0))/((-1.0+x)*x)

	solveZ   = xa .|> z

	if (solveZ[1] == Inf || isnan(solveZ[1]))
		solveZ[1] = 0.0
	end
	if (solveZ[lastindex(solveZ)] == Inf || isnan(solveZ[lastindex(solveZ)]))
		solveZ[lastindex(solveZ)] = 0.0
	end

	return 1.0/(NN2+0.0).*solveZ
end

function cumulativeSfs(sfsTemp)

	out    = Array{Float64}(undef, length(sfsTemp)+1)
	out[1] = sum(sfsTemp)

	for i in 2:length(sfsTemp)+1
		app = out[i-1]-sfsTemp[i-1]
		if app > 0.0
			out[i] = app
		else
			out[i] = (0.0)
		end
	end
	return out
end