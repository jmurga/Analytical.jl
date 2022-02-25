################################
######    Polymorphism    ######
################################
# Expected number of polymorphism above frequency x from the standard diffusion theory
# f(x) = ∫(s) θs* (1/x(1-x)) * ( (ℯ^(Ns)*(1-ℯ^(-4Ns(1-x))) ) / (ℯ^(4Ns-1)))
# Convolote with the binomial to obtain the downsampled sfs
# E[P(x)] = ∑_x=x*→x=1 fB(x*)
############Neutral#############

"""

	DiscSFSNeutDown()

Expected rate of neutral allele frequency reduce by backgrou	nd selection. The spectrum depends on the number of individual []

```math
\\mathbb{E}[Ps_{(x)}] = \\sum{x^{*}=x}{x^{*}=1}f_{B}(x)
```

# Return:
 - `Array{Float64}`: expected rate of neutral alleles frequencies.
"""
function DiscSFSNeutDown(param::parameters,binom::SparseMatrixCSC{Float64,Int64})

	NN2 = convert(Int64,ceil(param.NN*param.B))
	
	# Allocating variables
	neutral_sfs(i::Int64) = 1.0/(i)

	x = collect(0:NN2)
	solved_neutral_sfs = neutral_sfs.(x)
	replace!(solved_neutral_sfs, Inf => 0.0)

	# subsetDict = get(param.bn,param.B,1)
	# subsetDict = binom
	out::Array{Float64,1} = param.B * (param.θ_coding) * 0.25 * (binom*solved_neutral_sfs)

	return out
end

############Positive############
# Variable gamma in function changed to s to avoid problem with exported SpecialFunctions.gamma
"""

	DiscSFSSelPosDown(s,p)

Expected rate of positive selected allele frequency reduce by background selection. The spectrum depends on the number of individuals.

# Arguments
 - `s::Int64`: selection strength.
 - `p::Float64`: positive selected alleles probabilty.
# Return:
 - `Array{Float64}`: expected positive selected alleles frequencies.
"""
function DiscSFSSelPosDown(param::parameters,s::Int64,p::Float64,binom::SparseMatrixCSC{Float64,Int64})

	if p == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		# Solving sfs
		red_plus = Φ(param,s)
		NN2      = convert(Int64,ceil(param.NN*param.B))
		xa1      = collect(0:NN2)
		xa2      = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		s_corrected = s*param.B

		s_exp1 = exp(s_corrected*2)
		s_exp2 = exp(s_corrected*-2)

		positiveSfs(i::Float64,g1::Float64=s_exp1,g2::Float64=s_exp2,p::Float64=p) = Float64(p*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i))))

		# Original
		# p*0.5*(ℯ^(2*s_corrected)*(1-ℯ^(-2.0*s_corrected*(1.0-i)))/((ℯ^(2*s_corrected)-1.0)*i*(1.0-i)))

		# Allocating outputs
		solved_positive_sfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2))
		replace!(solved_positive_sfs, NaN => 0.0)

		# subsetDict = get(param.bn,param.B,1)
		# out               = (param.θ_coding)*red_plus*0.75*(subsetDict*solved_positive_sfs)
		out::Array{Float64,1} = param.θ_coding * red_plus * 0.75 * (binom*solved_positive_sfs)
		#=out::Array{Float64,1} = param.θᵣ .* red_plus .* 0.75 .* (binom*solved_positive_sfs)=#

	end

	return out
end

function DiscSFSSelPosDownArb(param::parameters,s::Int64,p::Float64,binom::SparseMatrixCSC{Float64,Int64})

	if p == 0.0
		out = zeros(Float64,param.nn + 1)
		out = out[2:end-1]
	else

		red_plus = Φ(param,s)

		# Solving sfs
		NN2  = convert(Int64,ceil(param.NN*param.B))
		xa1  = collect(0:NN2)
		xa2  = xa1/(NN2)

		# Solving float precision performance using exponential rule. Only one BigFloat estimation.
		s_corrected = s*param.B
		s_exp1::Quadmath.Float128 = exp(Quadmath.Float128(s_corrected*2))
		s_exp2::Quadmath.Float128 = exp(Quadmath.Float128(s_corrected*-2))

		positiveSfs(i::Float64,g1::Quadmath.Float128=s_exp1,g2::Quadmath.Float128=s_exp2,p::Float64=p) = Float64(p*0.5*(g1*(1- g2^(1.0-i))/((g1-1.0)*i*(1.0-i))))
		# Allocating outputs
		solved_positive_sfs::Array{Float64,1} = (1.0/(NN2)) * (positiveSfs.(xa2))
		replace!(solved_positive_sfs, NaN => 0.0)
		out::Array{Float64,1} = param.θ_coding * red_plus * 0.75 * (binom*solved_positive_sfs)
		#=out::Array{Float64,1} = param.θᵣ .* red_plus .* 0.75 .* (binom*solved_positive_sfs)=#

	end

	return out
end


######Slightly deleterious######
"""

	DiscSFSSelNegDown(param,p)

Expected rate of positive selected allele frequency reduce by background selection. Spectrum drawn on a gamma DFE. It depends on the number of individuals.

# Arguments
 - `p::Float64`: positive selected alleles probabilty.
# Return:
 - `Array{Float64}`: expected negative selected alleles frequencies.
"""
function DiscSFSSelNegDown(param::parameters,p::Float64,binom::SparseMatrixCSC{Float64,Int64})

	# subsetDict = get(param.bn,param.B,1)
	solved_negative = DiscSFSSelNeg(param,p)
	out = param.B .* (param.θ_coding) .* 0.75 .* (binom*solved_negative)

	return out
end

function DiscSFSSelNeg(param::parameters,p::Float64)

	beta     = param.scale/(1.0*param.B)
	NN2      = convert(Int64, ceil(param.NN*param.B))
	xa       = collect(0:NN2)/NN2

	solve_z   = similar(xa)

	z(x::Float64,p::Float64=p) = (1.0-p)*(2.0^-param.shape)*(beta^param.shape)*(-zeta(param.shape,x+beta/2.0) + zeta(param.shape,(2+beta)/2.0))/((-1.0+x)*x)

	solve_z   = z.(xa)

	if (solve_z[1] == Inf || isnan(solve_z[1]))
		solve_z[1] = 0.0
	end
	if (solve_z[lastindex(solve_z)] == Inf || isnan(solve_z[lastindex(solve_z)]))
		solve_z[lastindex(solve_z)] = 0.0
	end

	return 1.0/(NN2+0.0).*solve_z
end


"""
	cumulative_sfs(sfs_tmp)

Changing SFS considering all values above a frequency *x*. The original asymptotic-MK approach takes Pn(x) and Ps(x) as the number of polymorphic sites at frequency *x* rather than above *x*, but this approach scales poorly as sample size increases. We define the polymorphic spectrum as stated above since these quantities trivially have the same asymptote but are less affected by changing sample size.
"""
function cumulative_sfs(sfs_tmp::Array,freqs::Bool=true)

	out      = Array{Float64}(undef, size(sfs_tmp,1),size(sfs_tmp,2))

	if freqs
		idx = 2
	else
		idx = 1
	end

	out[1,idx:end] = sum(sfs_tmp[:,idx:end],dims=1)

	@simd for i in 2:(size(sfs_tmp)[1])

		#=app = view(out,i-1,:) .- view(sfs_tmp,i-1,:)=#
		app = out[i-1,idx:end] .- sfs_tmp[i-1,idx:end]

		if sum(app) > 0.0
			out[i,idx:end] = app
		else
			out[i,idx:end] = zeros(length(app))
		end
	end

	if freqs
		out[:,1] = sfs_tmp[:,1]
	end
	
	return out
end

"""
	reduce_sfs(sfs_tmp,bins)

Function to bin the SFS into a sample of N individuals.
"""
function reduce_sfs(sfs_tmp::Array,bins::Int64)

	n   = (bins*2) - 1 
	f   = sfs_tmp[:,1]
	sfs = sfs_tmp[:,2:end]
	
	b    = collect(1/n:1/n:1)
	inds = searchsortedfirst.(Ref(b), f)
	out  = zeros((n,size(sfs_tmp,2)))
	out[:,1] = unique(inds)
	sfs_grouped = hcat(inds,sfs)
	
	for i in unique(inds)
		out[out[:,1] .== i,2:end] = sum(sfs_grouped[sfs_grouped[:,1] .== i,2:end],dims=1)
	end
	out[:,1] = b

	return(out)
end
