using Distributions
using Statistics
using GZip
using CSV
using StatsBase

function readData(file)
	df = CSV.read(file)

	resampling = vcat([10,25,50,75],collect(100:100:1000),collect(1000:500:4000))
	tmp = zeros(length(resampling),2)
	for i in 1:length(resampling)
		idx = StatsBase.sample(axes(df, 1), resampling[i]; replace = true, ordered = true)
		tmp[i,:] = sum(convert(Array, df[idx,:]), dims = 1)
	end

	out = vcat(convert(Matrix,unique(df)),tmp,sum(convert(Array, df), dims = 1))
	return out
end


function meanQ(x)
	m = mean(x)
	q = Statistics.quantile(x[:,1],[0.05,0.95])
	return append!(q,m)
end

function ABCreg(;data::String, prior::String, nparams::Int64, nsummaries::Int64, outputPath::String, outputPrefix::String,tolerance::Float64, regressionMode::String,regPath="/home/jmurga/ABCreg/src/reg")

	reg = `$regPath -p $prior -d $data -P $nparams -S $nsummaries -b $outputPath/$outputPrefix -$regressionMode -t $tolerance`
	run(reg)

	files = filter(x -> occursin(outputPrefix,x), readdir(outputPath))
	# openFiles(f) = GZip.open(DelimitedFiles.readdlm,outputPath*"/"*f)
	openFiles(f) = convert(Matrix,CSV.read(GZip.open(outputPath*"/"*f)))

	estimates = Array{Float64}(undef,length(files),3)
	estimates = files .|> openFiles .|> meanQ
	results = reduce(hcat,estimates) |> transpose
	# results = Dict(zip(files, estimates))

	rm.(outputPath.*"/".*files)

	return results
end


# using Roots
# using HDF5
# module InferTools

# @with_kw mutable struct abcParams
# 	ABCreg::String         = ""
# end

# end #end module

# abcParameters = abcParams()
#
# function changeAbcParameters(;ABCreg="/home/jmurga/ABCreg/src/reg")
#
# 	abcParameters.ABCreg = ABCreg
#
# end


# function alphaF(d::Int64,d0::Int64,p::Int64,p0::Int64):
# 	if  d==0 or p0 == 0:
# 		return NaN
# 	al = 1.0 - (d0/d)*(p/p0)
# 	return al
# end
#
# function get_new_subs(num,alfac,befac,B)
# 	al = abcParameters.al
# 	be = abcParameters.be
# 	N = abcParameters.N
#
# 	z1 = SpecialFunctions.zeta(al,1+be/(2.0*B))
# 	z2 = SpecialFunctions.zeta(al,0.5*(2.0-1.0/N +be/B))
#
# 	denom = (z1-z2)
#
# 	z1n = SpecialFunctions.zeta(al*alfac,1+be*befac/(2.0*B))
# 	z2n = SpecialFunctions.zeta(al*alfac,0.5*(2.0-1.0/N +be*befac/B))
#
# 	numerator = (2.0^(-al*alfac+al)) * (B^(-al*alfac+al)) * (be^-al)* ((be*befac)^(al*alfac)) * (z1n-z2n)
#
# 	return convert(Int64,round(num*numerator/denom))
# end
#
# function bias_resample(arr,alfac,befac,al,be)
#
# 	gammaValues = (-1*arr)*abcParameters.N*2
#
# 	scale0 = pdf.(Distributions.Gamma(abcParameters.al,1/abcParameters.be),gammaValues)
# 	scale1 = pdf.(Distributions.Gamma(abcParameters.al*alfac,1/(abcParameters.be*befac)),gammaValues)
#
# 	weight = scale1./scale0
# 	weightret = sum(weight)/length(weight)
# 	weight = weight./sum(weight)
#
# 	samples = Distributions.sample(1:length(weight), Weights(weight),length(weight))
# 	return (samples,weightret)
# end
#
# function simulate()
# 	# now load all the B values for polymorphic sites, fixed sites
# 	Bh = readlines(Bfile)
#
# 	# Bdata is pN, pS, dN, dS
# 	id = collect(4:20)
# 	values = zeros(length(id),4)
# 	Bdata = Dict(id[i] => values[i,:] for i=1:length(id))
#
# 	for line in Bh
# 		data = split(strip(line))
# 		if parse(Float64,data[lastindex(data)])/1000. < 0.175
# 			continue
# 		end
# 		i = convert(Int64,ceil(20.0*parse(Float64,data[lastindex(data)])/1000.0))
#
# 		Bdata[i][1] += parse(Int64,data[2])
# 		Bdata[i][2] += parse(Int64,data[4])
# 		Bdata[i][3] += parse(Int64,data[6])
# 		Bdata[i][4] += parse(Int64,data[7])
# 	end
#
# 	# Loop: first sample theta_weak, theta_strong
# 	# For each B value in B_values, sample a polymorphic site with same B.
# 	# If site is NS, sample from NS frequency spectrum.
# 	sum_stats = [1,2,5,10,20,50,100,200,500,1000] .- 1
# end
#
# function installABCreg(path::String)
# 	clone = `git clone https://github.com/molpopgen/ABCreg.git $path/ABCreg`
# 	make = `make -C $path/ABCreg/src`
# 	run(clone)
# 	run(make)
#
# 	abcParameters.ABCreg = path * "/ABCreg/src/reg"
# end
#
