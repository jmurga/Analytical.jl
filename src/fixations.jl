################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutral fixation reduced by background selection

"""

	fixNeut()

Expected neutral fixations rate reduce by a background selection value. It takes into accoun the ammount probability of being synonymous. Testing build

```math
\\mathbb{E}[D_{s}] = (1 - p_{-} - p_{+}) B \\frac{1}{2N}
```

# Returns
	Expected rate of neutral fixations: Float64

"""
function fixNeut()
	return 0.255*(1.0/(adap.B*adap.NN))
end

# Negative fixations
"""

	fixNegB(ppos)

Expected fixation rate from negative DFE.

```math
\\mathbb{E}[D_{n-}] =  p_{-}\\left(2^-\\alpha\\beta^\\alpha\\left(-\\zeta[\\alpha,\\frac{2+\\beta}{2}] + \\zeta[\\alpha,1/2(2-\\frac{1}{N+\\beta})]\\right)\\right)
```

# Arguments
	ppos::Float64: selection coefficient

# Returns
	Expected rate of fixations from negative DFE: Float64

"""
function fixNegB(ppos::Float64)
	return 0.745*(1-ppos)*(2^(-adap.al))*(adap.B^(-adap.al))*(adap.be^adap.al)*(-SpecialFunctions.zeta(adap.al,1.0+adap.be/(2.0*adap.B))+SpecialFunctions.zeta(adap.al,0.5*(2-1.0/(adap.N*adap.B)+adap.be/adap.B)))
end

# Positive fixations
"""
	pFix()

Expected positive fixation rate
```math
\\mathbb{E}[D_{n+}] =  p_{+}  B (1 - 2^{(-2s)})
```

# Arguments
	ppos::Float64: selection coefficient

# Returns
	Rate of neutral fixation: Float64

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

Expected positive fixations rate reduced due to the impact of background selection and linkage. The formulas presented at the approach have been subjected to several previous theoretical works ([Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365), [Hudson et al., 1995](https://www.genetics.org/content/141/4/1605), (Nordborg et al. 1995)(https://doi.org/10.1017/S0016672300033619), [Barton NH., 1995](https://www.genetics.org/content/140/2/821)).

Then, the probabilty of fixation of positively selected alleles is reduced by a factor ```math \\phi``` across all deleterious linked sites. 


```math
\\phi(t,s) = \\euler^{[\\frac{-2\\mu}{}]}
```
```math
\\Phi = \\prod_{1}^{L} \\phi(t,s)
```
```math
\\mathbb{E}[D_{n+}]' =  \\Phi \\mathbb{E}[D_{n+}]
```
	
# Arguments
	ppos::Float64: selection coefficient
	
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

	# return 0.745 * ppos * phiReduction() * pFix(gamma)
	return 0.745 * ppos * ℯ^(-2.0*S*μ*(Ψ0-Ψ1)*CC^2/r^2) * pFix(gamma)
end
 