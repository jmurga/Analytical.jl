"""
	parse_sfs(;data,output,sfs_columns,div_columns)

Function to parse polymorphism and divergence by subset of genes. The input data is based on supplementary material described at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). Please be sure the file is tabulated.

| GeneId | Pn | DAF seppareted by commas | Ps | DAF separated by commas | Dn | Ds |
|--------|----|--------------------------|----|-------------------------|----|----|
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |

# Arguments
 - `data`: String or Array of strings containing files names with full path.
 - `output::String`: path to save file. Containing one file per input file.
 - `sfs_columns::Array{Int64,1}`: non-synonymous and synonymous daf columns. Please introduce first the non-synonymous number.
 - `div_columns::Array{Int64,1}`: non-synonymous and synonymous divergence columns. Please introduce first the non-synonymous number.

# Returns
 - `Array{Float64,1}`: α values
 - `Array{Float64,2}`: Site Frequency Spectrum
 - `Array{Float64,1}`: Synonymous and non-synonymous divergence counts
 - 
"""
function parse_sfs(;sample_size::Int64,data::String,gene_list::Union{Nothing,Any}=nothing,sfs_columns::Array{Int64,1}=[3,5],div_columns::Array{Int64,1}=[6,7],bins::Union{Nothing,Int64}=nothing,isolines::Bool=false)

	g(x) = parse.(Float64,x[2:end-1])
	
	if isolines
		s = sample_size
	else
		s = (sample_size*2)
	end
	
	freq = OrderedDict(round.(collect(1:(s-1))/s,digits=4) .=> 0)

	df   = CSV.read(data,header=false,delim='\t',DataFrame)

	if(!isnothing(gene_list))
		df =  vcat([ df[df[:,1] .==i,:]  for i in gene_list]...);
	end

	#=if(!isnothing(B))
		df = df[df[:,end] .== B,:]
		println(nrow(df))
		tmp  = split.(df[:,sfs_columns], ",")
	else
	end=#
	
	tmp  = split.(df[:,sfs_columns], ",")

	pn   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> countmap))
	ps   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> countmap))

	# Dn, Ds, Pn, Ps, sfs
	Dn           = sum(df[:,div_columns[1]])
	Ds           = sum(df[:,div_columns[2]])
	Pn           = sum(values(pn))
	Ps           = sum(values(ps))
	sfsPn        = reduce(vcat,values(merge(+,freq,pn)))
	sfsPs        = reduce(vcat,values(merge(+,freq,ps)))

	if(!isnothing(bins))
        sfsPn = reduce_sfs(hcat(collect(1:(s-1)),sfsPn),bins)[:,2]
        sfsPs = reduce_sfs(hcat(collect(1:(s-1)),sfsPs),bins)[:,2]

        sfs   = reduce_sfs(hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals),bins)
        scumu = cumulative_sfs(sfs)
	else
        sfs   = hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals)
        scumu = cumulative_sfs(sfs)
	end

    α    = round.(1 .- (Ds/Dn .*  scumu[:,2] ./scumu[:,3]),digits=5)

	return (α,sfs,[Dn,Ds])
end

"""
	ABCreg(analysis_folder, replicas, S, tol, abcreg)

Performing ABC inference using ABCreg. Please, be sure your analysis_folder contain the files alphas.txt and summaries.txt produced by Analytical.summaryStatsFromRates()

# Arguments
 - `analysis_folder::String` : Folder containing the observed data and summary estatistics. It will be used to output the posterior distributions
 - `S::Int64` : Number of summary stastitics to perform the inference.
 - `tol::Float64` : Tolerance value. It define the number of accepted value at ABC inference
 - `abcreg::String` : Path to ABCreg binary

# Output
Files containing posterior distributions from ABCreg

"""
function ABCreg(;analysis_folder::String,S::Int64,tol::Float64,abcreg::String)
	
	#=# List alphas and summstat files
    aFile   = analysis_folder * "/alphas.txt"
    sumFile = analysis_folder * "/summstat.txt"

	# Creating output names
	out = analysis_folder * "/out"

	r(a,s,o,abcreg=abcreg,S=S,tol=tol) = run(`$abcreg -d $a -p $s -P 5 -S $S -t $tol -b $o`)

	r(aFile,sumFile,out);=#

	# List alphas and summstat files
	aFile     = filter(x -> occursin("alphas",x), readdir(analysis_folder,join=true));
	sumFile   = filter(x -> occursin("summstat",x), readdir(analysis_folder,join=true));

	# Creating output names
	out = analysis_folder .* "/out_" .* string.(1:size(aFile,1))

	r(a,s,o,abcreg=abcreg,S=S,tol=tol) = run(`$abcreg -d $a -p $s -P 5 -S $S -t $tol -b $o`)

	progress_pmap(r,aFile,sumFile,out);	
end

"""
	Bootstrap data following polyDFE manual
"""
function bootstrap_data(sfs_files::Array{Float64,2},divergence_files::Array{Float64,2},replicas::Int64,output_folder::String)
	
	# Open Data
	sfs        = Array(CSV.read(sfs_files,DataFrame))
	divergence = fill(Array(CSV.read(divergence_files,DataFrame)),replicas)
	scumu      = fill(cumulative_sfs(sfs[:,2:end]),replicas)

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

function open_sfs_div(x::Array{String,1},y::Array{String,1},dac::Array{Int64,1},replicas::Int64,bootstrap::Bool)

	sfs = Array.(CSV.read.(x,DataFrame,header=false))
	divergence = Array.(CSV.read.(y,DataFrame,header=false))

	if bootstrap
		sfs = repeat(sfs,replicas)
		divergence = repeat(divergence,replicas)
		pr(x) = hcat(x[:,1],pois_rand.(x[:,2:end]))
		sfs[2:end] .= pr.(sfs[2:end])
	end

	scumu = cumulative_sfs.(sfs)
	f(x,d=dac) = sum(x[:,2:3],dims=2)[d]
	s = f.(scumu)

	d = [[sum(divergence[i][1:2])] for i in eachindex(divergence)]
	al(a,b,c=dac) = @. round(1 - (b[2]/b[1] * a[:,2]/a[:,3])[c],digits=5)
	α = permutedims.(al.(scumu,divergence))

	return(s,d,α,map(x -> x[:,2:end],scumu))
end


function install_r()

		@eval using Pkg
		const ENV["R_HOME"]="*"

		@eval Pkg.add("Conda")

		@eval using Conda
		@eval Conda.add("r-base",channel="conda-forge")
		@eval Conda.add(["r-locfit","r-ggplot2","r-data.table","r-r.utils"],channel="conda-forge")

		@eval Pkg.add("RCall")
end
"""
	Function to download and source plotMap function. We do not include at part of the module to avoid external dependecies. Once the function is execute properly you will have a function called *plot_map which used R to 
		estimate and plot the Maximum A Posterior following ABCreg example. It uses locfit and ggplot2 libraries.
"""
#=function source_plot_map_R(;script::String)=#
#=function plot_map(;analysis_folder::String,weak::Bool=true,title::String="Posteriors")

	try
		@eval R"""library(ggplot2);library(locfit);library(data.table)"""

		out = filter(x -> occursin("post",x), readdir(analysis_folder,join=true))
		out          = filter(x -> !occursin(".1.",x),out)

		open(x)      = Array(CSV.read(x,DataFrame,header=false))

		# Control outlier inference. 2Nes non negative values
		flt(x)       = x[x[:,4] .> 0,:]
		posteriors   = flt.(open.(out))

		@eval R"""getmap <- function(df){
				temp = as.data.frame(df)
				d <-locfit(~temp[,1],temp);
				map<-temp[,1][which.max(predict(d,newdata=temp))]
			}"""

		getmap(x)    = rcopy(R"""suppressWarnings(matrix(apply($x,2,getmap),nrow=1))""")
		
		if !weak
			posteriors = [posteriors[i][:,3:end] for i in eachindex(posteriors)]
			tmp          = getmap.(posteriors)
			maxp         = DataFrame(vcat(tmp...),[:a,:gamNeg,:shape])
			al           = maxp[:,1:1]
			gam          = maxp[:,2:end]
		else
			tmp          = getmap.(posteriors)
			maxp         = DataFrame(vcat(tmp...),[:aw,:as,:a,:gamNeg,:shape])
			al           = maxp[:,1:3]
			gam          = maxp[:,4:end]
		end

		@eval R"""al = as.data.table($al)
			lbls = if(ncol(al) > 1){c(expression(paste('Posterior ',alpha[w])), expression(paste('Posterior ',alpha[s])),expression(paste('Posterior ',alpha)))}else{c(expression(paste('Posterior ',alpha)))}
			clrs = if(ncol(al) > 1){c('#30504f', '#e2bd9a', '#ab2710')}else{c('#ab2710')}


			dal = suppressWarnings(data.table::melt(al))
			pal = suppressWarnings(ggplot(dal) + geom_density(aes(x=value,fill=variable),alpha=0.75) + scale_fill_manual('Posterior distribution',values=clrs ,labels=lbls) + theme_bw() + ggtitle($title) + xlab(expression(alpha)) + ylab(""))
			suppressWarnings(ggsave(pal,filename=paste0($analysis_folder,'/map.png'),dpi=600))
			"""
		CSV.write(analysis_folder * "/map.tsv",maxp,delim='\t',header=true)

		#=RCall.endEmbeddedR();=#

		return(maxp)
	catch
		println("Please install R, ggplot2, data.table and locfit in your system before execute this function")
		install_r()
	end
end
=#