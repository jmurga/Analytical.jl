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
 - `Array{Float64,1}`: α values
 - `Array{Float64,2}`: Site Frequency Spectrum
 - `Array{Float64,1}`: Synonymous and non-synonymous divergence counts
 - 
"""
function parseSfs(;sampleSize::Int64,data::String,geneList::Union{Nothing,Array{String,1}}=nothing,sfsColumns::Array{Int64,1}=[3,5],divColumns::Array{Int64,1}=[6,7],bins::Union{Nothing,Int64}=nothing)

	g(x) = parse.(Float64,x[2:end-1])
	
	s = (sampleSize*2)
	freq = OrderedDict(round.(collect(1:(s-1))/s,digits=4) .=> 0)

	df   = CSV.read(data,header=false,delim='\t',DataFrame)

	if(!isnothing(geneList))
		df = df[∈(geneList).(df[:,1]), :]
	end

	#=if(!isnothing(B))
		df = df[df[:,end] .== B,:]
		println(nrow(df))
		tmp  = split.(df[:,sfsColumns], ",")
	else
	end=#
	
	tmp  = split.(df[:,sfsColumns], ",")

	pn   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> countmap))
	ps   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> countmap))

	# Dn, Ds, Pn, Ps, sfs
	Dn           = sum(df[:,divColumns[1]])
	Ds           = sum(df[:,divColumns[2]])
	Pn           = sum(values(pn))
	Ps           = sum(values(ps))
	sfsPn        = reduce(vcat,values(merge(+,freq,pn)))
	sfsPs        = reduce(vcat,values(merge(+,freq,ps)))

	if(!isnothing(bins))
        sfsPn = reduceSfs(hcat(collect(1:(s-1)),sfsPn),bins)[:,2]
        sfsPs = reduceSfs(hcat(collect(1:(s-1)),sfsPs),bins)[:,2]

        sfs   = reduceSfs(hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals),bins)
	else
        sfs   = hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals)
        scumu = cumulativeSfs(sfs)
	end

    α    = round.(1 .- (Ds/Dn .*  scumu[:,2] ./scumu[:,3]),digits=5)

	return (α,sfs,[Dn,Ds])
end

"""
	parseBinnedSfs(;data,output,sfsColumns,divColumns)

Function to parse polymorphism and divergence by subset of genes. The input data is based iMKT [Murga-Moreno et. al 2019](https://doi.org/10.1038/s41559-019-0890-6). Please be sure the file is tabulated.

| GeneId | Pn | DAF0f separated by ; | Ps | DAF0F separated by ; | Dn | Ds |
|--------|----|--------------------------|----|-------------------------|----|----|
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |

# Arguments
 - `data`: String or Array of strings containing files names with full path.

# Returns
 - `Array{Int64,1}`: α values
 - `Array{Int64,1}`: Cumulative SFS
 - `Array{Int64,2}`: Total divergence counts
 - `Array{Int64,1}`: Synonymous and non-synonymous divergence counts
 - 
"""
function parseBinnedSfs(;data::String,population::String,geneList::Union{Nothing,Array{String,1}}=nothing,cumulative::Bool=false)

	
    df   = CSV.read(data,header=true,delim='\t',DataFrame)

    tmp  = df[df.pop .== population, : ]
    tmp  = tmp[(tmp.pi .!=0) .&(tmp.p0 .!=0),:]
    tmp  = tmp[(tmp.di .!=0) .&(tmp.d0 .!=0),:]
	
	if(!isnothing(geneList))
		tmp = tmp[∈(geneList).(tmp.id), :]
	end

	ds   = sum(tmp.d0)
	dn   = sum(tmp.di)

	ms   = sum(tmp.m0)
	mn   = sum(tmp.mi)

	s(x) = parse.(Int,split(x,";"))

	sfsPs = sum(reduce(hcat,s.(tmp.daf4f)),dims=2)
	sfsPn = sum(reduce(hcat,s.(tmp.daf0f)),dims=2)

	if cumulative
		f = collect(1:size(sfsPs,1)) ./ (size(sfsPs,1)+1)
		sfs = cumulativeSfs(hcat(f,sfsPn,sfsPs))
	else
		f = collect(1:size(sfsPs,1)) ./ (size(sfsPs,1)+1)
		sfs = hcat(f,sfsPn,sfsPs)
	end

	α = @. 1 - (ds/dn * sfs[:,2]/sfs[:,3])
	return(sfs,[dn,ds],[mn,ms],α)
end

"""
	ABCreg(analysis, replicas, P, S, tol, workers, abcreg, parallel)

Could be parallelize if GNU parallel is available in your system
"""
function ABCreg(;analysis::String,replicas::Int64,P::Int64,S::Int64,tol::Float64,workers::Int64,abcreg::String,parallel::Bool=false)
	
	# List alphas and summstat files
	aFile = @. analysis * "/alpha_" * string(1:replicas) * ".tsv"
	sumFile = @. analysis * "/summstat_" * string(1:replicas) * ".tsv"

	# Creating output names
	out = @. analysis * "/out_" * string(1:replicas)

	# Check if GNU parallel installed
	if parallel
		run(`parallel -j$workers --link $abcreg -d "{1}" -p "{2}" -P $P -S $S -t $tol -b "{3}" ::: $aFile ::: $sumFile ::: $out`)
	else
		r(a,s,o) = run(`$abcreg -d $a -p $s -P $P -S $S -t $tol -b $o`)
		r.(aFile,sumFile,out)
	end	
end


"""
	Bootstrap data following polyDFE manual
"""
function bootstrapData(sFile::Array{Float64,2},dFile::Array{Float64,2},replicas::Int64,outputFolder::String)
	
	# Open Data
	sfs        = Array(CSV.read(sFile,DataFrame))
	divergence = fill(Array(CSV.read(dFile,DataFrame)),replicas)
	scumu      = fill(cumulativeSfs(sfs[:,2:end]),replicas)

	# Bootstraping
	b(x)       = pois_rand.(x)
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

	sfs = Array.(CSV.read.(x,DataFrame,header=false))
	divergence = Array.(CSV.read.(y,DataFrame,header=false))

	if bootstrap
		sfs = repeat(sfs,replicas)
		divergence = repeat(divergence,replicas)
		pr(x) = hcat(x[:,1],pois_rand.(x[:,2:end]))
		sfs = pr.(sfs)
	end


	scumu = Analytical.cumulativeSfs.(sfs)
	f(x,d=dac) = sum(x[:,2:3],dims=2)[d]
	s = f.(scumu)

	d = [[sum(divergence[i])] for i in eachindex(divergence)]
	al(a,b,c=dac) = @. round(1 - (b[2]/b[1] * a[:,2]/a[:,3])[c],digits=5)
	α = permutedims.(al.(scumu,divergence))
	return(s,d,α)	
end

"""
	Function to download and source plotMap function. We do not include at part of the module to avoid external dependecies. Once the function is execute properly you will have a function called *plotMap which used R to 
		estimate and plot the Maximum A Posterior following ABCreg example. It uses locfit and ggplot2 libraries.
"""
function sourcePlotMapR(;script::String)

	download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/plotMapR.jl",script)

	try
		include(script)
	catch
		"Please be sure you installed RCall, R and abc library in your system. Check our documentation if you have any installation problem: https://jmurga.github.io/Analytical.jl/dev/"
	end
end