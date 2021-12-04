################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutral fixation corrected by background selection
"""

 - fixNeut()

Expected neutral fixations rate reduce by B value.

```math
\\mathbb{E}[D_{s}] = (1 - p_{-} - p_{+}) B \\frac{1}{2N}
```
# Returns
 - `Float64`: expected rate of neutral fixations.

"""
function fixNeut(param::parameters)
	# Synonymous probabilty * (fixation probabilty corrected by BGS value)
	out::Float64 = 0.25*(1.0/(param.B*param.NN))
	return out
end

# Negative fixations corrected by background selection
"""

	fixNegB(ppos)

Expected fixation rate from negative DFE.

```math
\\mathbb{E}[D_{n-}] =  p_{-}\\left(2^-\\alpha\\beta^\\alpha\\left(-\\zeta[\\alpha,\\frac{2+\\beta}{2}] + \\zeta[\\alpha,1/2(2-\\frac{1}{N+\\beta})]\\right)\\right)
```

# Arguments
 - `ppos::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of fixations from negative DFE.

"""
function fixNegB(param::parameters,ppos::Float64)
	# Non-synonymous proportion * negative alleles probability * fixation probability from gamma distribution
	out::Float64 = 0.75*(1-ppos)*(2^(-param.al))*(param.B^(-param.al))*(param.be^param.al)*(-zeta(param.al,1.0+param.be/(2.0*param.B))+zeta(param.al,0.5*(2-1.0/(param.N*param.B)+param.be/param.B)))
	return out
end

# Positive fixations
"""
	pFix()

Expected positive fixation rate.

```math
\\mathbb{E}[D_{n+}] =  p_{+} \\cdot B \\cdot (1 - e^{(-2s)})
```

# Arguments
 - `ppos::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of positive fixation.

"""
function pFix(param::parameters,gam::Int64)

	# Fixation probability
	s::Float64    = gam/(param.NN)
	pfix::Float64 = (1.0-ℯ^(-2.0*s))/(1.0-ℯ^(-2.0*gam))

	# Correcting pFix for large s following Uricchio et al. 2014
	if s >= 0.1
		pfix = ℯ^(-(1.0+s))
		lim = 0
		while(lim < 200)
			pfix = ℯ^((1.0+s)*(pfix-1.0))
			lim  = lim + 1
		end
		pfix = 1 -pfix
	end

	return pfix
end

# Positive fixations after apply Φ. reduction of positive fixations due deleterious linkage given a value B of background selection
"""

	fixPosSim(gamma,ppos)

Expected positive fixations rate reduced due to the impact of background selection and linkage. The probabilty of fixation of positively selected alleles is reduced by a factor Φ across all deleterious linked sites [`Analytical.phiReduction`](@ref).

```math
\\mathbb{E}[D_{n+}] =  \\Phi \\cdot \\mathbb{E}[D_{n+}]
```

# Arguments
 - `ppos::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of positive fixations under background selection

"""
function fixPosSim(param::parameters,gamma::Int64,ppos::Float64)

	redPlus = phiReduction(param,gamma)

	# Non-synonymous * positive alleles probability * B reduction * fixation probility
	out::Float64 = 0.75*ppos*redPlus*pFix(param,gamma)
	return out
end
