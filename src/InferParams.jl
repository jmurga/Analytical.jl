module InferTools

using Parameters
# using PyCall
# using SpecialFunctions
using Distributions
using StatsBase
# using Roots
# using HDF5

@with_kw mutable struct abcParams

	Bfile::String          = ""
	al::Float64            = 0.184
	be::Float64            = 0.000402
	N::Int64               = 500
	datadir::String        = "~/ABC/output"
	datadirtempr::String   = "~/ABC/output_temp"
	prefr::String          = "test"
	nsims::Int64           = 1000
	ABCsumsdir::String     = "ABC_sums_neg_test"
	alfac::Float64         = 1.0
	befac::Float64         = 1.0
	testsim::Bool          = false
	reg::Bool              = true
	q::Int64               = -1
	print_als::Bool        = false
	gerp_control::Bool     = false
	Brange::Array{Float64} = Array{Float64}(undef,10)
	task_id::Int64         = 1
end

abcParameters = abcParams()


Bfile  = test.Bfile
nsims = test.nsims
if 'SLURM_ARRAY_TASK_ID' in os.environ:
	test.task_id = int(os.environ['SLURM_ARRAY_TASK_ID'])

function alphaF(d::Int64,d0::Int64,p::Int64,p0::Int64):
	if  d==0 or p0 == 0:
		return NaN
	al = 1.0 - (d0/d)*(p/p0)
	return al
end

function get_new_subs(num,alfac,befac,B)
	al = abcParameters.al
	be = abcParameters.be
	N = abcParameters.N

	z1 = SpecialFunctions.zeta(al,1+be/(2.0*B))
	z2 = SpecialFunctions.zeta(al,0.5*(2.0-1.0/N +be/B))

	denom = (z1-z2)

	z1n = SpecialFunctions.zeta(al*alfac,1+be*befac/(2.0*B))
	z2n = SpecialFunctions.zeta(al*alfac,0.5*(2.0-1.0/N +be*befac/B))

	numerator = (2.0^(-al*alfac+al)) * (B^(-al*alfac+al)) * (be^-al)* ((be*befac)^(al*alfac)) * (z1n-z2n)

	return convert(Int64,round(num*numerator/denom))
end

function bias_resample(arr,alfac,befac,al,be)

	gammaValues = (-1*arr)*abcParameters.N*2

	scale0 = pdf.(Distributions.Gamma(abcParameters.al,1/abcParameters.be),gammaValues)
	scale1 = pdf.(Distributions.Gamma(abcParameters.al*alfac,1/(abcParameters.be*befac)),gammaValues)

	weight = scale1./scale0
	weightret = sum(weight)/length(weight)
	weight = weight./sum(weight)

	samples = Distributions.sample(1:length(weight), Weights(weight),length(weight))
	return (samples,weightret)
end

function simulate()
	# now load all the B values for polymorphic sites, fixed sites
	Bh = readlines(Bfile)

	# Bdata is pN, pS, dN, dS
	id = collect(4:20)
	values = zeros(length(id),4)
	Bdata = Dict(id[i] => values[i,:] for i=1:length(id))

	for line in Bh
		data = split(strip(line))
		if parse(Float64,data[lastindex(data)])/1000. < 0.175
			continue
		end
		i = convert(Int64,ceil(20.0*parse(Float64,data[lastindex(data)])/1000.0))

		Bdata[i][1] += parse(Int64,data[2])
		Bdata[i][2] += parse(Int64,data[4])
		Bdata[i][3] += parse(Int64,data[6])
		Bdata[i][4] += parse(Int64,data[7])
	end

	# Loop: first sample theta_weak, theta_strong
	# For each B value in B_values, sample a polymorphic site with same B.
	# If site is NS, sample from NS frequency spectrum.
	sum_stats = [1,2,5,10,20,50,100,200,500,1000] .- 1
end

end #end module
