"""
	parseSfs(;data,output,sfsColumns,divColumns)

Function to parse polymorphism and divergence by subset of genes. The input data is based on supplementary material described at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). Please be sure the file is tabulated.

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
function parseSfs(;sample::Int64,data::String,sfsColumns::Array{Int64,1}=[3,5],divColumns::Array{Int64,1}=[6,7],dac::Array{Int64,1},B::Union{Nothing,Float64}=nothing,bins::Union{Nothing,Int64}=nothing)

	g(x) = parse.(Float64,x[2:end-1])
	
	s = (sample*2)
	freq = OrderedDict(round.(collect(1:(s-1))/s,digits=4) .=> 0)

	df   = CSV.read(data,header=false,delim='\t',DataFrame)
	#=df   = filter([:Column2, :Column4] => (x, y) -> x > 0 || y > 0 , df)=#
	

	if(!isnothing(B))
		df = df[df[:,end] .== B,:]
		println(nrow(df))
		tmp  = split.(df[:,sfsColumns], ",")
	else
		tmp  = split.(df[:,sfsColumns], ",")
	end

	pn   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> StatsBase.countmap))
	ps   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> StatsBase.countmap))

	# Dn, Ds, Pn, Ps, sfs
	Dn           = sum(df[:,divColumns[1]])
	Ds           = sum(df[:,divColumns[2]])
	Pn           = sum(values(pn))
	Ps           = sum(values(ps))
	sfsPn        = cumulativeSfs(reduce(vcat,values(merge(+,freq,pn))))
	sfsPs        = (cumulativeSfs(reduce(vcat,values(merge(+,freq,ps)))))

	if(!isnothing(bins))
		sfsPn = reduceSfs(hcat(collect(1:(s-1)),sfsPn),40)[:,2]
		sfsPs = reduceSfs(hcat(collect(1:(s-1)),sfsPs),40)[:,2]

		sfs = reduceSfs(hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals),bins)
	else
		sfs = hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals)
	end
	α            = round.(1 .- Ds/Dn .*  sfsPn ./sfsPs,digits=5)[dac]

	D = [Dn+Ds]
	cSfs = sfsPn+sfsPs

	return (α,cSfs,D,sfs,[Dn,Ds])
end
"""
	ABCreg(data, prior, nparams, nsummaries, outputPath, outputPrefix,tolerance, regressionMode,regPath)
 -
"""
function ABCreg(;analysis::String,replicas::Int64,P::Int64,S::Int64,tol::Float64,workers::Int64,parallel::Bool)
	
	reg = chomp(read(`which reg`,String))

	aFile = @. analysis * "/alpha_" * string(1:replicas) * ".tsv"
	sumFile = @. analysis * "/summstat_" * string(1:replicas) * ".tsv"
	out = @. analysis * "/out_" * string(1:replicas)

	if parallel
		run(`parallel -j$workers --link $reg -d "{1}" -p "{2}" -P $P -S $S -t $tol -L -b "{3}" ::: $aFile ::: $sumFile ::: $out`)
	else
		r(a,s,o) = run(`reg -d $a -p $s -P $P -S $S -t $tol -L -b $o`)
		r.(aFile,bFile,out)
	end	
end


"""
	Bootstrap data
"""
function bootstrapData(sFile::Array{Float64,2},dFile::Array{Float64,2},replicas::Int64,outputFolder::String)
	
	# Open Data
	sfs        = Array(CSV.read(sFile,DataFrame))
	divergence = fill(Array(CSV.read(dFile,DataFrame)),replicas)
	scumu      = fill(cumulativeSfs(sfs[:,2:end]),replicas)

	# Bootstraping
	b(x)       = PoissonRandom.pois_rand.(x)
	bootstrap  = b.(scumu)

	# Output
	outSfs = @. output * "sfs" * string(1:replicas) * ".tsv"
	outDiv = @. output * "div" * string(1:replicas) * ".tsv"
	for i  = 1:replicas
		tmp = hcat(sfs[:,1],bootstrap[i])
		CSV.write(out,DataFrame(tmp),header=false)
		CSV.write(out,DataFrame(divergence[i]),header=false)
	end
end

function openSfsDiv(x::Array{String,1},y::Array{String,1},dac::Array{Int64,1},replicas::Int64,bootstrap::Bool)

	sfs = Array.(CSV.read.(x,DataFrame))
	divergence = Array.(CSV.read.(y,DataFrame))

	if bootstrap
		sfs = repeat(sfs,replicas)
		divergence = repeat(divergence,replicas)
		pr(x) = hcat(x[:,1],PoissonRandom.pois_rand.(x[:,2:end]))
		sfs = pr.(sfs)
	end


	scumu = cumulativeSfs.(sfs)
	f(x,d=dac) = sum(x[:,2:3],dims=2)[d]
	s = f.(scumu)

	d = [[sum(divergence[i])] for i in eachindex(divergence)]
	al(a,b,c=dac) = @. round(1 - (b[2]/b[1] * a[:,2]/a[:,3])[c],digits=5)
	α = al.(scumu,divergence)
	return(s,d,α)	
end

"""
	Estimating and plotting MAP using locfit and ggplot2 in R. It assume your folder contains the posterior estimated through ABCreg
"""
function sourcePlotMap(;outputFile::String)

	download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/plotMapR.jl",outputFile)

	include(outputFile)

end