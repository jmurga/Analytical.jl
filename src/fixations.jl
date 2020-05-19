################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutral fixation reduced by background selection

"""

 - fixNeut()

Expected neutral fixations rate reduce by a background selection value.

```math
\\mathbb{E}[D_{s}] = (1 - p_{-} - p_{+}) B \\frac{1}{2N}
```
# Returns
 - `Float64`: expected rate of neutral fixations.

"""
function fixNeut()
			# SYNONYMOUS * NEGATIVE PROBABILITY * FIXATION PROBABILITY FROM GAMMA DISTRIBUTION
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
 - `ppos::Float64`: selection coefficient.

# Returns
 - `Float64`: expected rate of fixations from negative DFE.

"""
function fixNegB(ppos::Float64)
		# NON-SYNONYMOUS * NEGATIVE PROBABILITY * FIXATION PROBABILITY FROM GAMMA DISTRIBUTION
	return 0.745*(1-ppos)*(2^(-adap.al))*(adap.B^(-adap.al))*(adap.be^adap.al)*(-SpecialFunctions.zeta(adap.al,1.0+adap.be/(2.0*adap.B))+SpecialFunctions.zeta(adap.al,0.5*(2-1.0/(adap.N*adap.B)+adap.be/adap.B)))
end

# Positive fixations
"""
	pFix()

Expected positive fixation rate.

```math
\\mathbb{E}[D_{n+}] =  p_{+} \\cdot B \\cdot (1 - e^{(-2s)})
```

# Arguments
 - `ppos::Float64`: selection coefficient.

# Returns
 - `Float64`: expected rate of positive fixation.

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

Expected positive fixations rate reduced due to the impact of background selection and linkage. The probabilty of fixation of positively selected alleles is reduced by a factor Φ across all deleterious linked sites [`Analytical.phiReduction`](@ref).

```math
\\mathbb{E}[D_{n+}] =  \\Phi \\cdot \\mathbb{E}[D_{n+}]
```
	
# Arguments
 - `ppos::Float64`: selection coefficient
	
# Returns
 - `Float64`: expected rate of positive fixations under background selection

"""
function fixPosSim(gamma::Int64,ppos::Float64)

	# S  = abs(adap.gam_neg/(1.0*adap.NN))
	# r  = adap.rho/(2.0*adap.NN)
	# μ  = adap.theta_f/(2.0*adap.NN)
	# s  = gamma/(adap.NN*1.0)

	# Ψ0 = SpecialFunctions.polygamma(1,(s+S)/r)
	# Ψ1 = SpecialFunctions.polygamma(1,(r+adap.Lf*r+s+S)/r)
	# CC = 1.0
	red_plus = phiReduction(gamma)

		# NON-SYNONYMOUS * POSITIVE PROBABILITY * BGS REDUCTION * FIXATION PROB
	return 0.745 * ppos * red_plus * pFix(gamma)
end

