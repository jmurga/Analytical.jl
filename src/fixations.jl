################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutral fixation reduced by background selection
"""
	fixNeut()

Expected neutral fixations rate reduce by a background selection value. It takes into accoun the ammount probability of being synonymous.

```math	
\\mathbb{E}\\left[D_{s}\\rigth] =  \\left(1 - p_{-} - p_{+}) \\cdot B \\cdot \\frac{1}{2N}\\rigth)
````
# Returns
	Expected rate of neutral fixations: Float64
"""
function fixNeut()
	return 0.255*(1.0/(adap.B*adap.NN))
end

# Negative fixations
"""

	fixNegB(ppos)

Expected fixation rate from negative DFE
```math
\\mathbb{E}\\left[D_{n-}\\right] =  p_{-} \\left(2^-\\alpha \\cdot \\beta^\\alpha \\cdot \\left(-\\zeta\\left[\\alpha,2+\\beta/2] + \\zeta\\left(\\left[\\alpha,1/2*\\left(\\left(2-\\fract{1}{N+\\beta}\\rigth)\\rigth]\\rigth)\\rigth)
```
# Arguments
	- ```ppos::Float64```: selection coefficient
Negative fixations.
# Returns
	Expected rate of fixations from negative DFE in a Float64
"""
function fixNegB(ppos::Float64)
	return 0.745*(1-ppos)*(2^(-adap.al))*(adap.B^(-adap.al))*(adap.be^adap.al)*(-SpecialFunctions.zeta(adap.al,1.0+adap.be/(2.0*adap.B))+SpecialFunctions.zeta(adap.al,0.5*(2-1.0/(adap.N*adap.B)+adap.be/adap.B)))
end

# Positive fixations
"""
	fixNegB(pFix)

Expected positive fixation rate
```math
\\mathbb{E}\\left[D_{s}\\rigth] =  \\left(1 - p_{-} - p_{+}) \\cdot B \\cdot \\frac{1}{2N}\\rigth)
```
# Arguments
	- ```ppos::Float64```: selection coefficient
# Returns
	Rate of neutral fixation in a Float64
"""
function pFix(gamma::Int64)

	s    = gamma/(adap.NN+0.0)
	pfix = (1.0-ℯ^(-2.0*s))/(1.0-ℯ^(-2.0*gamma))

	if s >= 0.1
		pfix = ℯ^(-(1.0+s))
		lim = 0
		while(lim < 200)
			pfix = ℯ^((1.0+s)*(pfix-1.0))
			lim +=1
		pfix = 1-pfix
		end
	end

	return pfix
end

# Positive fixations after apply Φ, reduction of positive fixations due deleterious linkage given a value B of background selection
"""

	fixPosSim(gamma,ppos)

Reduction of expected positive fixations rate due deleterious linkage given a value  of background selection. The fixation probability of positively selected alleles is reduced by a factor ```math \\phi``` across all deleterious linked sites

```math
	\\mathbb{E}\\left[D_{s}\\rigth] = p_{+}\\left(1-\\euler^-2s\\right)
	\\phi\\left(t,s\\right) = e^\\left(\\frac{-2\\mu}{t\\left(1+\\frac{rL}{t}+\\frac{2s}{t}\\right)}\\rigth)
	\\Phi = \\prod_{1}^{L} \\phi(t,s)
	\\mathbb{E}\\left[D_{s}'\\rigth] = \\mathbb{E}\\left[D_{s}\\rigth] \\cdot \\Phi
```

	
# Arguments
	- ```ppos::Float64```: selection coefficient
# Returns
	Expected rate of positive fixations under background selection: Float64
"""
function fixPosSim(gamma::Int64,ppos::Float64)

	S  = abs(adap.gam_neg/(1.0*adap.NN))
	r  = adap.rho/(2.0*adap.NN)
	μ  = adap.theta_f/(2.0*adap.NN)
	s  = gamma/(adap.NN*1.0)

	Ψ0 = SpecialFunctions.polygamma(1,(s+S)/r)
	Ψ1 = SpecialFunctions.polygamma(1,(r+adap.Lf*r+s+S)/r)
	CC = 1.0

	# return 0.745 * ppos * phiReduction * pFix(gamma)
	return 0.745 * ppos * ℯ^(-2.0*S*μ*(Ψ0-Ψ1)*CC^2/r^2) * pFix(gamma)
end
 