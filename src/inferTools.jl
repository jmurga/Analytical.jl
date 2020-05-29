"""
	parseSfs(;data,output,sfsColumns,divColumns)

Function to parse polymorphism and divergence by subset of genes. The input data is based on supplementary material described at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6).

| GeneId | Pn | DAF seppareted by commas | Ps | DAF separated by commas | Dn | Ds |
|--------|----|--------------------------|----|-------------------------|----|----|
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |

# Arguments
 - `data`: String or Array of strings containing files names with full path.
 - `output::String`: path to save file. Containing one file per input file.
 - `sfsColumns::Array{Int64,1}`: non-synonymous and synonymous daf columns. Please introduce first the non-synonymous number.
 - `divColumns::Array{Int64,1}`: non-synonymous and synonymous divergence columns. Please introduce first the non-synonymous number.
 
# Returns
 - `Array{Array{Int64,N} where N,1}`: Array of arrays containing the total polymorphic sites (1), total Site Frequency Spectrum (2) and total divergence (3). Each array contains one row/column per file.
 - File writed in `output`
"""
# function parseSfs(;data,output::String,sfsColumns::Array{Int64,1}=[3,5],polColumns::Array{Int64,1}=[2,4],divColumns::Array{Int64,1}=[6,7],bins::Int64)
	
# 	g(x) = parse.(Float64,x[2:end-1])

# 	function binned(tmpDf,b::Int64)
		
# 		tmp1   = round.(reduce(vcat,tmpDf[:,1] .|> g),digits=4)
# 		tmp2   = round.(reduce(vcat,tmpDf[:,2] .|> g),digits=4)
# 		tmp = vcat(tmp1,tmp2)
		
# 		tmpOrdered = DataFrame(hcat(tmp,cut(tmp,0:(1/b):1;extend=true))) |> sort
# 		binnedSfs = combine(nrow,groupby(tmpOrdered, :x2))

# 		return binnedSfs[:,2]
# 	end

# 	if(data isa String)

# 		P   = Array{Int64}(undef,1)
# 		D   = Array{Int64}(undef,1)
# 		sfs = Array{Float64}(undef, bins ,1)
# 		newData = Array{Float64}(undef, 1,24)

# 		df = read(data,header=false,delim=' ')
# 		pol  = convert(Matrix,split.(df[:,sfsColumns], ","))

# 		# Empirical data to analytical estimations
# 		sfs = binned(pol,bins) |> permutedims
# 		P = sum(sfs)
# 		D = convert(Matrix,df[:,divColumns]) |> sum
# 		newSfs = binned(pol,20) |> permutedims
# 		# Dn, Ds, Pn, Ps, sfs
# 		newData = [sum(df[:,divColumns[1]]) sum(df[:,divColumns[2]]) sum(df[:,polColumns[1]]) sum(df[:,polColumns[2]]) newSfs]
		
# 		write(output, DataFrame(newData), delim='\t',writeheader=false)
# 		return [P,sfs,D]
# 	else

# 		P   = Array{Float64}(undef,length(data))
# 		D   = Array{Float64}(undef,length(data))
# 		sfs = Array{Float64}(undef,length(data),bins)
# 		newData = Array{Int64}(undef,length(data),24)

# 		for i in 1:length(data)

# 			df = read(data[i],header=false,delim=' ')
# 			pol  = convert(Matrix,split.(df[:,sfsColumns], ","))

# 			# Empirical data to analytical estimations
# 			sfs[i,:] = binned(pol,bins) |> permutedims
# 			P[i] = sum(sfs[i,:])
# 			D[i] = convert(Matrix,df[:,divColumns]) |> sum
# 			newSfs = binned(pol,20) |> permutedims
# 			# Dn, Ds, Pn, Ps, sfs
# 			newData[1,:] = [sum(df[:,divColumns[1]]) sum(df[:,divColumns[2]]) sum(df[:,polColumns[1]]) sum(df[:,polColumns[2]]) newSfs]

# 		end

# 		write(output, DataFrame(newData),delim='\t',writeheader=false)
# 		return [P,permutedims(sfs),D]
# 	end
# end
function parseSfs(;data,output::String,sfsColumns::Array{Int64,1}=[3,5],divColumns::Array{Int64,1}=[6,7],bins::Int64)
	
	g(x) = parse.(Float64,x[2:end-1])
	
	freq = OrderedDict(round.(collect(1:adap.nn-1)/adap.nn,digits=4) .=> 0)

	if(data isa String)

		P   = Array{Int64}(undef,1)
		D   = Array{Int64}(undef,1)
		sfs = Array{Float64}(undef, adap.nn -1 ,1)
		newData = Array{Float64}(undef, 1,24)

		df = read(data,header=false,delim=' ')
		df = filter([:Column2, :Column4] => (x, y) -> x > 0 || y > 0 , df)
		tmp  = split.(df[:,sfsColumns], ",")
		pn   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> StatsBase.countmap))
		ps   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> StatsBase.countmap))


		# x = zeros(adap.nn)
		# y = zeros(adap.nn)
		# for i in 1:adap.nn
		# 	round.((i/adap.nn),digits=4) 
		# 	try
		# 		x[i] = pn[round.((i/adap.nn),digits=4)]
		# 		y[i] = ps[round.((i/adap.nn),digits=4)]
		# 	catch
		# 		x[i] = 0
		# 		y[i] = 0
		# 	end
		# end

		# Empirical data to analytical estimations
		tmpSfs   =  merge(+,pn,ps)
		sfs = reduce(vcat,values(merge(+,freq,tmpSfs)))
		P = sum(sfs)
		D = convert(Matrix,df[:,divColumns]) |> sum
		# Dn, Ds, Pn, Ps, sfs
		newData = [sum(df[:,divColumns[1]]) sum(df[:,divColumns[2]]) sum(values(pn)) sum(values(ps)) reduceSfs(sfs,bins)]
		
		write(output, DataFrame(newData), delim='\t',writeheader=false)
		return [P,sfs,D]

	else
		P   = Array{Int64}(undef,length(data))
		D   = Array{Int64}(undef,length(data))
		sfs = Array{Float64}(undef,length(data),adap.nn-1)
		newData = Array{Int64}(undef,length(data),4+bins)
		
		for i in 1:length(data)
			
			df = read(data[i],header=false,delim=' ')
			df = filter([:Column2, :Column4] => (x, y) -> x > 0 || y > 0 , df)
			tmp  = split.(df[:,sfsColumns], ",")
			pn   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> StatsBase.countmap))
			ps   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> StatsBase.countmap))

			# Empirical data to analytical estimations
			tmpSfs   =  merge(+,pn,ps)
			sfs[i,:] = reduce(vcat,values(merge(+,freq,tmpSfs)))
			P[i]  = sum(sfs[i,:][1])
			D[i] = convert(Matrix,df[:,divColumns]) |> sum

			newData[i,:] = [sum(df[:,divColumns[1]]) sum(df[:,divColumns[2]]) sum(values(pn)) sum(values(ps)) reduceSfs(sfs[i,:],bins)]

		end

		write(output, DataFrame(newData),delim='\t',writeheader=false)
		return [P,permutedims(sfs),D]

	end
end


"""
	ABCreg(data, prior, nparams, nsummaries, outputPath, outputPrefix,tolerance, regressionMode,regPath)

Julia to execute *ABCreg*. You should take into account that any other ABC could be used once the prior distributions are done. If you consider to use ABCreg please [cite the publication](https://doi.org/10.1186/1471-2156-10-35) and compile it in your system.

# Arguments
 - `data::String`: Observed data. Produced by parseSfs.
 - `output::String`: path to save file. ABCreg will produce one file per lines inside data.
 - `nparams::Int64`: number of parameters in prior distributions.
 - `nsummaries::Int64`: number of observed summary statistics.
 - `outputPath::String`: output path.
 - `outputPrefix::String`: output prefix name.
 - `tolerance::Float64`: tolerance for rejection acceptance in Euclidian distance.
 - `regressionMode::String`: transformation (T or L).
# Returns
 - `Array{Array{Int64,N} where N,1}`: Array of arrays containing the total polymorphic sites (1), total Site Frequency Spectrum (2) and total divergence (3). Each array contains one row/column per file.
 - File writed in `output`
 - 
"""
function ABCreg(;data::String, prior::String, nparams::Int64, nsummaries::Int64, outputPath::String, outputPrefix::String,tolerance::Float64, regressionMode::String,regPath="/home/jmurga/ABCreg/src/reg")

	reg = `$regPath -p $prior -d $data -P $nparams -S $nsummaries -b $outputPath$outputPrefix -$regressionMode -t $tolerance`

	openFiles(f) = convert(Matrix,read(open(f),header=false))

	run(reg)

	files = outputPath .* filter(x -> occursin(outputPrefix,x), readdir(outputPath))

	posteriors = files .|> openFiles
	estimates  = posteriors .|> meanQ
	
	return posteriors,estimates
end

"""
	meanQ(x,columns)

Function to retrieve mean and quantiles (95%) from posterior distributions.

# Arguments
 - `x::Array{Float64,2}`: posterior distribution.
 - `columns::Array{Int64,1}`: columns to process.
# Returns
 - `Array{Array{Float64,2},1}`: Array of array containing mean and quantiles by posterior distribution. Each array contains ```\$\\alpha_{S}\$```, ```\$\\alpha_{W}\$``` and ```\$\\alpha\$``` information by column.
"""
function meanQ(x::Array{Float64,2})

	m            = StatsBase.mean(x,dims=1)
	qt           = Array{Float64}(undef,size(x,2),2)
	for i in 1   : size(x,2)
		qt[i,: ] = StatsBase.quantile(x[:,i],[0.05,0.95])
	end

	return vcat(m,permutedims(qt))
end

function readData(file)
	df = read(file)

	resampling = vcat([10,25,50,75],collect(100:100:1000),collect(1000:500:4000))
	tmp = zeros(length(resampling),2)
	for i in 1:length(resampling)
		idx = StatsBase.sample(axes(df, 1), resampling[i]; replace = true, ordered = true)
		tmp[i,:] = sum(convert(Array, df[idx,:]), dims = 1)
	end

	out = vcat(convert(Matrix,unique(df)),tmp,sum(convert(Array, df), dims = 1))
	return out
end

function reduceSfs(sfsTemp,bins)

	freq  = collect(0:size(sfsTemp,1)-1)/size(sfsTemp,1)
	h1    = fit(Histogram,freq,0:(1/bins):1)
	xmap1 = StatsBase.binindex.(Ref(h1), freq)
	
	
	tmp = hcat(sfsTemp,xmap1)
	out = zeros(bins,size(sfsTemp,2))
	for i in unique(xmap1)
		out[i,:] = sum(tmp[tmp[:,end].==i,1:end-1],dims=1)
	end
		
	return (permutedims(out))
end