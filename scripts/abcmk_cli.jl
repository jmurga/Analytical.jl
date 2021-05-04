using Fire, Distributed

"Function to estimate rates based on prior values. Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function rates(;ne::Int64=1000, samples::Int64=500, gamNeg::String="-1000 -200", gL::String="5 10", gH::String="400 1000",dac::String="2,4,5,10,20,50,200,500,700",shape::Float64=0.184,rho::String="nothing",theta::String="nothing",solutions::Int64=100000,output::String="/home/jmurga/rates.jld2",nthreads::Int64=1)

	tmpNeg    = parse.(Int,split(gamNeg," "))
	tmpStrong = parse.(Int,split(gH," "))
	dac       = parse.(Int,split(dac,","))

	if (gL == "nothing")
		tmpWeak = nothing
	else
		tmpWeak = parse.(Int,split(gL," "))
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

	addprocs(nthreads)
	
	@eval @everywhere using Analytical
	@eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac,al=$shape)

	@eval convolutedSamples = Analytical.binomialDict()
	@eval Analytical.binomOp!($adap,$convolutedSamples.bn);
	@time @eval df = Analytical.rates(param = $adap,convolutedSamples=$convolutedSamples,gH=collect($tmpStrong[1]:$tmpStrong[2]),gL=$tmpWeak,gamNeg=collect($tmpNeg[1]:$tmpNeg[2]),iterations = $solutions,rho=$rho,theta=$theta,shape=$adap.al,output=$output);

	# remove the workers
	for i in Distributed.workers()
		rmprocs(i)
	end
end

"Function to parse polymorphic and divergence data from Uricchio et. al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/

	julia abcmk_cli.jl parseTgpData --analysisFolder /home/jmurga/test/abcmk/dna2/ --geneList /home/jmurga/test/abcmk/dnaVipsList.txt

"
@main function parseTgpData(;analysisFolder::String="<folder>",geneList::String="false",bins::String="false")
	
	@eval using Analytical, DataFrames, CSV, ProgressMeter

	run(`mkdir -p $analysisFolder`)

	tgpData = analysisFolder * "/mk_with_positions_BGS.txt"
	
	download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/mk_with_positions_BGS.txt",tgpData)

	# Check if bins or genelist are defined
	@eval if $geneList != "false"
		@eval gList = CSV.read($geneList,DataFrame,header=false)[:,1]
	else
		gList = nothing
	end

	@eval if $bins != "false"
		@eval binsSize = parse(Int,$bins)
	else
		binsSize = nothing
	end

	# Parsing TGP data
	@eval α,sfs, divergence = Analytical.parseSfs(sampleSize=661,data=$tgpData,geneList=$gList,bins=$binsSize)

	# Writting data to folder
	@eval sName = $analysisFolder * "/sfs.tsv"
	@eval dName = $analysisFolder * "/div.tsv"

	@eval CSV.write($sName,DataFrame($sfs),delim='\t',header=false)
	@eval CSV.write($dName,DataFrame($divergence'),delim='\t',header=false)
end

"Please provide an analysis folder containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function parseDgnData(;analysisFolder::String="<analysisFolder>",geneList::String="false",population::String="ZI")
	
	@eval using Analytical, DataFrames, CSV, ProgressMeter
	@eval run(`mkdir -p $analysisFolder`)

end


"Estimate summary statistics from analytical rates. You must provide a path containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function summaries(;analysisFolder::String="<folder>",rates::String="rates.jld2",ne::Int64=1000, samples::Int64=500,dac::String="2,4,5,10,20,50,200,500,700",summstatSize::Int64=100000,replicas::Int64=100,bootstrap::String="true",nthreads::Int64=1)
	
	#=ne=1000;samples=661;rates="/home/jmurga/test.jld2";summstatSize=100000;replicas=100;analysisFolder="/home/jmurga/test/abcmk/dna2/";bootstrap="true";dac="2,4,5,10,20,50,200,661,921"=#

	addprocs(nthreads)
	
	@eval @everywhere using Analytical, DataFrames, CSV, ProgressMeter
	@eval using JLD2
	
	@eval h5file    = jldopen($rates)

	@eval adap      = Analytical.parameters(N=$ne,n=$samples,dac =parse.(Int,split($dac,",")))

	@eval if $bootstrap == "true"
		@eval summstat  = Analytical.summaryStatsFromRates(param=$adap,rates=$h5file,analysisFolder=$analysisFolder,summstatSize=$summstatSize,replicas=$replicas,bootstrap=true)
	else
		@eval summstat  = Analytical.summaryStatsFromRates(param=$adap,rates=$h5file,analysisFolder=$analysisFolder,summstatSize=$summstatSize,replicas=$replicas,bootstrap=false)
	end
end

"ABCreg inference"
@main function abcInference(;analysisFolder::String="<folder>",replicas::Int64=100,P::Int64=5,S::Int64=9,tol::Float64=0.01,nthreads::Int64=1,ABCreg::String="/home/jmurga/ABCreg/src/reg",parallel::String="false")
	
	@eval aFile   = @. $analysisFolder * "/alpha_" * string(1:$replicas) * ".tsv"
	@eval sumFile = @. $analysisFolder * "/summstat_" * string(1:$replicas) * ".tsv"
	@eval out     = @. $analysisFolder * "/out_" * string(1:$replicas)

	if parallel  == "true"
		run(`parallel -j$nthreads --link $ABCreg -d "{1}" -p "{2}" -P 5 -S 9 -t $tol -b "{3}" ::: $aFile ::: $sumFile ::: $out`);
	else
		@eval r(a,s,o) = run(`reg -d $a -p $s -P 5 -S 9 -t 0.01 -b $o`)
		@eval r.($aFile,$sumFile,$out);
	end	
end

"Plot Maximum a posterior distribution"
@main function plotMap(;analysisFolder::String="<folder>",output::String="<folder>")

	try
		@eval using RCall, GZip, DataFrames, CSV
		@eval Analytical.sourcePlotMapR(script=$analysisFolder * "/script.jl")
		@eval plotMap(analysisFolder=$analysisFolder);
	catch
		println("Please install R, ggplot2 and abc in your system before execute this function")
	end
end

