using Fire, Distributed

"Function to estimate rates based on prior values. Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
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

	if scheduler == "local"
		@eval addprocs($nthreads)
	elseif scheduler == "slurm"
		@eval using ClusterManagers
		@eval addprocs_slurm($nthreads)
	elseif scheduler == "htcondor"
		@eval using ClusterManagers
		@eval addprocs_htc($nthreads)
	end
	
	@eval @everywhere using Analytical
	@eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac,al=$shape)

	@eval convolutedSamples = Analytical.binomialDict()
	@eval Analytical.binomOp!($adap,$convolutedSamples.bn);
	@eval df = Analytical.rates(param = $adap,convolutedSamples=$convolutedSamples,gH=collect($tmpStrong[1]:$tmpStrong[2]),gL=$tmpWeak,gamNeg=collect($tmpNeg[1]:$tmpNeg[2]),iterations = $solutions,rho=$rho,theta=$theta,shape=$adap.al,output=$output);
	for i in workers()
		rmprocs(i)
	end
end

"Function to parse polymorphic and divergence data from Uricchio et. al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/

	julia abcmk_cli.jl parseTgpData --analysisFolder /home/jmurga/test/abcmk/dna2/ --geneList /home/jmurga/test/abcmk/dnaVipsList.txt
"
@main function parseData(;analysisFolder::String="<folder>",dataset::String="tgp",geneList::String="false",bins::String="false")
	
	@eval using Analytical, DataFrames, CSV

	run(`mkdir -p $analysisFolder`)

	data = analysisFolder * "/" * dataset * ".txt"
	
	download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/"* dataset * ".txt",data)

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
	@eval α,sfs, divergence = Analytical.parseSfs(sampleSize=661,data=$data,geneList=$gList,bins=$binsSize)

	# Writting data to folder
	@eval sName = $analysisFolder * "/sfs.tsv"
	@eval dName = $analysisFolder * "/div.tsv"

	@eval CSV.write($sName,DataFrame($sfs,:auto),delim='\t',header=false)
	@eval CSV.write($dName,DataFrame($divergence',:auto),delim='\t',header=false)
end

"Estimate summary statistics from analytical rates. You must provide a path containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function summaries(;analysisFolder::String="<folder>",rates::String="rates.jld2",ne::Int64=1000, samples::Int64=500,dac::String="2,4,5,10,20,50,200,500,700",summstatSize::Int64=100000,replicas::Int64=100,bootstrap::String="true",scheduler::String="local",nthreads::Int64=1)
	

	if scheduler == "local"
		@eval addprocs($nthreads)
	elseif scheduler == "slurm"
		@eval using ClusterManagers
		@eval addprocs_slurm($nthreads)
	elseif scheduler == "htcondor"
		@eval using ClusterManagers
		@eval addprocs_htc($nthreads)
	end
	
	@eval @everywhere using Analytical, ParallelUtilities
	@eval using JLD2, DataFrames, CSV, ProgressMeter
	
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
        @eval r(a,s,o,r=$ABCreg,nP=$P,nS=$S,t=$tol) = run(`$r -d $a -p $s -P $nP -S $nS -t $t -b $o`)
        @eval r.($aFile,$sumFile,$out);
	end	
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

