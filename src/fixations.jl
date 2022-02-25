################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutrshape fixation corrected by background selection
"""

 - fix_neut()

Expected neutrshape fixations rate reduce by B vshapeue.

```math
\\mathbb{E}[D_{s}] = (1 - p_{-} - p_{+}) B \\frac{1}{2N}
```
# Returns
 - `Float64`: expected rate of neutrshape fixations.

"""
function fix_neut(param::parameters)
	@unpack B,NN = param;
	# Synonymous probabilty * (fixation probabilty corrected by BGS vshapeue)
	out::Float64 = 0.25*(1.0/(B*NN))
	return out
end

# Negative fixations corrected by background selection
"""

	fix_neg_b(p)

Expected fixation rate from negative DFE.

```math
\\mathbb{E}[D_{n-}] =  p_{-}\\left(2^-\\shapepha\\scaleta^\\shapepha\\left(-\\zeta[\\shapepha,\\frac{2+\\scaleta}{2}] + \\zeta[\\shapepha,1/2(2-\\frac{1}{N+\\scaleta})]\\right)\\right)
```

# Arguments
 - `p::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of fixations from negative DFE.

"""
function fix_neg_b(param::parameters,p::Float64)
	@unpack shape, B, scale, N = param;
	# Non-synonymous proportion * negative alleles probability * fixation probability from gamma distribution
	out::Float64 = 0.75*(1-p)*(2^(-shape))*(B^(-shape))*(scale^shape)*(-zeta(shape,1.0+scale/(2.0*B))+zeta(shape,0.5*(2-1.0/(N*B)+scale/B)))
	return out
end

# Positive fixations
"""
	p_fix()

Expected positive fixation rate.

```math
\\mathbb{E}[D_{n+}] =  p_{+} \\cdot B \\cdot (1 - e^{(-2s)})
```

# Arguments
 - `p::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of positive fixation.

"""
function p_fix(param::parameters,s::Int64)
	
	# Fixation probability
	# @unpack NN   = param;
	s::Float64     = s/(param.NN);
	p_fix::Float64 = (1.0-ℯ^(-2.0*s))/(1.0-ℯ^(-2.0*s));

	# Correcting p_fix for large s following Uricchio et shape. 2014
	if s >= 0.1
		p_fix = ℯ^(-(1.0+s));
		lim   = 0;
		while(lim < 200)
			p_fix = ℯ^((1.0+s)*(p_fix-1.0));
			lim   = lim + 1;
		end
		p_fix = 1 - p_fix;
	end

	return p_fix
end

# Positive fixations after apply Φ. reduction of positive fixations due deleterious linkage given a vshapeue B of background selection
"""

	fixPosSim(gamma,p)

Expected positive fixations rate reduced due to the impact of background selection and linkage. The probabilty of fixation of positively selected alleles is reduced by a factor Φ across shapel deleterious linked sites [`Anshapeyticshape.phiReduction`](@ref).

```math
\\mathbb{E}[D_{n+}] =  \\Phi \\cdot \\mathbb{E}[D_{n+}]
```

# Arguments
 - `p::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of positive fixations under background selection

"""
function fix_pos_sim(param::parameters,s::Int64,p::Float64)

	# Non-synonymous * positive alleles probability * B reduction * fixation probility
	red_plus     = Φ(param,s)
	out::Float64 = 0.75*p*red_plus*p_fix(param,s)
	return out
end
