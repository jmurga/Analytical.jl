using Fire, Distributed

"
	julia abcmk_cli.jl rates --samples 661 --gamNeg -1000,-200 --gL 1,10 --gH 400,1000 --rho 0.001 --theta 0.001 --solutions 1000000 --output rates.jld2 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000  --nthreads 7

Function to solve fixation and polymorphic rates analitically. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.

If rho and/or theta are set to nothing, the function will input random values given the range 0.0005:0.0005:0.01. Otherwise you can fix the values.
If gL is set to nothing, the function will not account the role of the weakly selected alleles in the estimation.

The function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.

Please check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/Analytical.jl/dev/cli/ to check model
"
@main function rates(;ne::Int64=1000, samples::Int64=500, gamNeg::String="-1000,-200", gL::String="5,10", gH::String="400,1000",dac::String="2,4,5,10,20,50,200,500,700",shape::Float64=0.184,rho::String="nothing",theta::String="nothing",solutions::Int64=100000,output::String="/home/jmurga/rates.jld2",scheduler::String="local",nthreads::Int64=1)

	tmpNeg    = parse.(Int,split(gamNeg,","))
	tmpStrong = parse.(Int,split(gH,","))
	dac       = parse.(Int,split(dac,","))

	if (gL == "nothing")
		tmpWeak = nothing
	else
		tmpWeak = parse.(Int,split(gL,","))
		tmpWeak = collect(tmpWeak[1]:tmpWeak[2])
	end


	if (rho == "nothing")
		rho = nothing
	else
		rho = parse(Float64,rho)
	end

	if (theta == "nothing")
		theta = nothing
	else
		theta = parse(Float64,theta)
	end

	if scheduler == "slurm"
		@eval using ClusterManagers
		@eval addprocs_slurm($nthreads)
	elseif scheduler == "htcondor"
		@eval using ClusterManagers
		@eval addprocs_htc($nthreads)
	else
		@eval addprocs($nthreads)
	end
	
	@eval @everywhere using Analytical, ParallelUtilities
	@eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac,al=$shape)

	@eval convolutedSamples = Analytical.binomialDict()
	@eval Analytical.binomOp!($adap,$convolutedSamples.bn);
	@eval Analytical.rates(param = $adap,convolutedSamples=$convolutedSamples,gH=collect($tmpStrong[1]:$tmpStrong[2]),gL=$tmpWeak,gamNeg=collect($tmpNeg[1]:$tmpNeg[2]),iterations = $solutions,rho=$rho,theta=$theta,shape=$adap.al,output=$output,scheduler=$scheduler);

	for i in workers()
		rmprocs(i)
	end
end


"
	julia abcmk_cli.jl joinRates --analysisFolder --output
"
@main function joinRates(;analysisFolder::String="<folder>",samples::Int64=661,output::String="<name>")

	@eval using JLD2, CSV ,DataFrames, OrderedCollections

	searchdir(path) = filter(x->occursin("jld2",x), readdir(path))
	files = analysisFolder .* "/" .* searchdir(analysisFolder)

	j = jldopen.(files)
	
	models = []
	n = []
	s = []
	dsdn = []
	dac = j[1]["1000/" * string(samples)]["dac"]


	for i in eachindex(files)
		tmp = j[i]["1000/" * string(samples)]

		push!(models,tmp["models"])
		push!(n,tmp["neut"])
		push!(s,tmp["sel"])
		push!(dsdn,tmp["dsdn"])

	end

	neut = OrderedDict{Int,Array}()
	sel = OrderedDict{Int,Array}()
	o(x,k) = x[k]
	for x in keys(n[1])
		neut[x] = vcat(o.(n,x)...)
		sel[x] = vcat(o.(s,x)...)
	end
	models = vcat(models...)
	dsdn = vcat(dsdn...)

	# Writting HDF5 file
	JLD2.jldopen(output, "a+") do file
		file[string(1000)* "/" * string(samples) * "/models"] = models;
		file[string(1000)* "/" * string(samples) * "/neut"]   = neut;
		file[string(1000)* "/" * string(samples) * "/sel"]    = sel;
		file[string(1000)* "/" * string(samples) * "/dsdn"]   = dsdn;
		file[string(1000)* "/" * string(samples) * "/dac"]    = dac;
	end;

	rm.(files)
end


"
	julia abcmk_cli.jl parseData --analysisFolder analysis/ --geneList analysis/dnaVipsList.txt

Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. 

The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics.	

Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli/
"
@main function parseData(;analysisFolder::String="<folder>",dataset::String="tgp",geneList::String="false",bins::String="false")
	
	@eval using Analytical, DataFrames, CSV

	run(`mkdir -p $analysisFolder`)

	dataset = lowercase(dataset)
	data    = analysisFolder * "/" * dataset * ".txt"
	
	download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/"* dataset * ".txt",data)

	# Check if bins or genelist are defined
	@eval if $geneList != "false"
		@eval gList = CSV.read($geneList,DataFrame,header=false) |> Array
	else
		gList = nothing
	end

	@eval if $bins != "false"
		@eval binsSize = parse(Int,$bins)
	else
		binsSize = nothing
	end

	# Parsing TGP data
	if dataset == "tgp"
		@eval α,sfs, divergence = Analytical.parseSfs(sampleSize=661,data=$data,geneList=$gList,bins=$binsSize)
	elseif occursin("zi",dataset)
		@eval α,sfs, divergence = Analytical.parseSfs(sampleSize=154,data=$data,geneList=$gList,bins=$binsSize,isolines=true)
	elseif occursin("ral",dataset)
		@eval α,sfs, divergence = Analytical.parseSfs(sampleSize=160,data=$data,geneList=$gList,bins=$binsSize,isolines=true)
	end
	# Writting data to folder
	@eval sName = $analysisFolder * "/sfs.tsv"
	@eval dName = $analysisFolder * "/div.tsv"

	@eval CSV.write($sName,DataFrame($sfs,:auto),delim='\t',header=false)
	@eval CSV.write($dName,DataFrame($divergence',:auto),delim='\t',header=false)
end

"
	julia abcmk_cli.jl summaries --analysisFolder analysis/ --rates analysis/rates.jld2 --samples 661 --dac 2,4,5,10,20,50,200,661,925 --summstatSize 1000000 --nthreads 7

Estimate summary statistics from analytical rates. You must provide a path containing the parsed SFS and divergence file.

The function returns files containing bootstrapped datasets (alphas.txt) and summary statistics (summstat.txt)

Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli"
@main function summaries(;analysisFolder::String="<folder>",rates::String="rates.jld2",ne::Int64=1000, samples::Int64=500,dac::String="2,4,5,10,20,50,200,500,700",summstatSize::Int64=100000,replicas::Int64=100,bootstrap::String="true")
	
	@eval  using Analytical, JLD2, DataFrames, CSV, ProgressMeter
	
	@eval h5file    = jldopen($rates)

	@eval adap      = Analytical.parameters(N=$ne,n=$samples,dac =parse.(Int,split($dac,",")))

	@eval if $bootstrap == "true"
		@eval summstat  = Analytical.summaryStatsFromRates(param=$adap,rates=$h5file,analysisFolder=$analysisFolder,summstatSize=$summstatSize,replicas=$replicas,bootstrap=true)
	else
		@eval summstat  = Analytical.summaryStatsFromRates(param=$adap,rates=$h5file,analysisFolder=$analysisFolder,summstatSize=$summstatSize,replicas=$replicas,bootstrap=false)
	end
end

"ABCreg inference.

The function returns posterior distributions from ABC inference. Each posterior file contains information about alpha_w, alpha_s, alpha, gamNeg and shape parameter. The number of posterior distributions will depend on the number of bootstrap replicas.

Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli
"
@main function abcInference(;analysisFolder::String="<folder>",S::Int64=9,tol::Float64=0.01,ABCreg::String="/home/jmurga/ABCreg/src/reg")
	
	@eval using Analytical
	@eval Analytical.ABCreg(analysisFolder=$analysisFolder,S=$S,tol=$tol,abcreg=$ABCreg)

end

#="Plot Maximum a posterior distribution"
@main function plotMap(;analysisFolder::String="<folder>")
	try
		@eval using Analytical, RCall, GZip, DataFrames, CSV
		
		@eval Analytical.sourcePlotMapR(script=$analysisFolder)
		@eval Analytical.plotMap(analysisFolder=$analysisFolder)
		@eval RCall.endEmbeddedR()
	catch
		println("Please install R, ggplot2, data.table and locfit in your system before execute this function")
	end
end=#

