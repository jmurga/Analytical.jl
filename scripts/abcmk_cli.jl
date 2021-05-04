using Fire, Distributed

"Function to estimate rates"
@main function rates(;ne::Int64=1000, samples::Int64=500, gamNeg::String="-1000 -200", gL::String="5 10", gH::String="400 1000",dac::String="2,4,5,10,20,50,200,500,700",shape::Float64=0.184,rho::String="nothing",theta::String="nothing",solutions::Int64=1000000,output::String="/home/jmurga/rates.jld2",workers::Int64=1)

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

	addprocs(parse(Int,workers))
	

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

"Please provide an analysis folder containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function parseTgpData(;analysisFolder::String="<folder>",geneList::String="false",bins::String="false")
	
	@eval using Analytical, DataFrames, CSV, ProgressMeter

	@eval run(`mkdir -p $analysisFolder`)

	tgpData = analysisFolder * "/mk_with_positions_BGS.txt"
	
	@eval download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/mk_with_positions_BGS.txt",$tgpData)

	# Check if bins or genelist are defined
	@eval if geneList != "false"
		gList = CSV.read(geneList,DataFrame,header=false)[:,1]
	else
		gList = nothing
	end

	@eval if bins != "false"
		binsSize = parse(Int,$bins)
	else
		binsSize = nothing
	end

	# Parsing TGP data
	@eval α,sfs, divergence = parseSfs(sampleSize=661,data=$tgpData,geneList=$gList,bins=$binsSize)

	# Writting data to folder
	sName = analysisFolder * "/sfs.tsv"
	dName = analysisFolder * "/div.tsv"

	@eval CSV.write($sName,DataFrame($sfs),delim='\t',header=false)
	@eval CSV.write($dName,DataFrame($divergence'),delim='\t',header=false)

end

"Please provide an analysis folder containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function parseDgnData(;analysisFolder::String,geneList::String="false",population::String="ZI")
	
	@eval using Analytical, DataFrames, CSV, ProgressMeter
	@eval run(`mkdir -p $analysisFolder`)

end



"Summary statistics from rates. Please provide an analysis folder containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function summaries(;analysisFolder::String="<folder>",rates::String="rates.jld2",ne::Int64=1000, samples::Int64=500,dac::String="2,4,5,10,20,50,200,500,700",summstatSize::Int64=100000,replicas::Int64=100,bootstrap::String="true",nthreads::Int64=1)
	
	#=ne=1000;samples=661;rates="/home/jmurga/test.jld2";summstatSize=100000;replicas=100;analysisFolder="/home/jmurga/test/abcmk/dna/";bootstrap="true";dac="2,4,5,10,20,50,200,661,921"=#

	addprocs(nthreads)
	
	@eval @everywhere using Analytical, DataFrames, CSV, ProgressMeter
	@eval using JLD2
	
	@eval h5file    = jldopen($rates)

	@eval adap      = Analytical.parameters(N=$ne,n=$samples,dac =parse.(Int,split($dac,",")))

	@eval if $bootstrap == "true"
		@eval summstat  = Analytical.summaryStatsFromRates(;param=$adap,rates=$h5file,analysisFolder=$analysisFolder,summstatSize=$summstatSize,replicas=$replicas,bootstrap=true)
	else
		@eval summstat  = Analytical.summaryStatsFromRates(;param=$adap,rates=$h5file,analysisFolder=$analysis,summstatSize=$summstatSize,replicas=$replicas,bootstrap=false)
	end
end



"ABCreg inference"
@main function abcInference(;analysis::String="<folder>",replicas::Int64=100,P::Int64=5,S::Int64=9,tol::Float64=0.01,workers::Int64=1,parallel::String="true")
	
	reg           = chomp(read(`which reg`,String))
	@eval aFile   = @. $analysis * "/alpha_" * string(1:$replicas) * ".tsv"
	@eval sumFile = @. $analysis * "/summstat_" * string(1:$replicas) * ".tsv"
	@eval out     = @. $analysis * "/out_" * string(1:$replicas)
	if parallel  == "true"
		run(`parallel -j$workers --link $reg -d "{1}" -p "{2}" -P 5 -S 9 -t 0.01 -L -b "{3}" ::: $aFile ::: $sumFile ::: $out`)
	else
		@eval r(a,s,o) = run(`reg -d $a -p $s -P 5 -S 9 -t 0.01 -L -b $o`)
		@eval r.($aFile,$sumFile,$out)
	end	
end

"Plot Maximum a posterior distribution"
@main function plotMap(;analysis::String="<folder>",output::String="<folder>")

	try
		@eval using RCall, GZip, DataFrames, CSV
		@eval R"""library(ggplot2);library(abc)"""
		@eval R"""getmap <- function(df){
					temp = as.data.frame(df)
					d <-locfit(~temp[,1],temp);
					map<-temp[,1][which.max(predict(d,newdata=temp))]}"""

		out              = filter(x -> occursin("post",x), readdir(analysis,join=true))
		out              = filter(x -> !occursin(".1.",x),out)

		open(x)          = Array(CSV.read(GZip.open(x),DataFrame,header=false))
		@eval posteriors = open.($out)
		
		@eval getmap(x)  = rcopy(R"""matrix(apply(x,2,getmap),nrow=1)""")
		@eval tmp        = getmap.($posteriors)
		@eval maxp       = DataFrame($tmp,[:aw,:as,:a,:gamNeg,:shape])
	
		al               = maxp[:,1:3]
		gam              = maxp[:,4:end]

		@eval @rput(al)
		@eval @rput(output)
		@eval p = R"""al = al
			dal = melt(al)
			pal = ggplot(dal) + geom_density(aes(x=value,fill=variable),alpha=0.5) + scale_fill_manual(values=c('#30504f', '#e2bd9a', '#ab2710'))
			ggsave(pal,filename=paste0(output,'/posteriorAlphas.svg'))"""
			
		@eval CSV.write(output * "/mapDistribution.tsv",maxp,header=true,delim='\t')
	catch
		println("Please install R, ggplot2 and abc in your system before execute this function")
	end
end

