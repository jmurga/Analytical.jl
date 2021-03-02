using Fire, Distributed

"Function to estimate rates"
@main function rates(;ne::Int64=1000, samples::Int64=500, gamNeg::String="-1000 -200", gL::String="5 10;nothing", gH::String="400 1000",dac::String="2,4,5,10,20,50,200,500,700",shape::Float64=0.184,rho::Float64=0.001,theta::Float64=0.001,solutions::Int64=1000000,output::String="/home/jmurga/rates.jld2",workers::Int64=1,cluster::String="local")

	tmpNeg    = parse.(Int,split(gamNeg," "))
	tmpStrong = parse.(Int,split(gH," "))
	dac       = parse.(Int,split(dac,","))

	if (tmpWeak == "nothing")
		gL = nothing
	else
		gL = parse.(Int,split(gL," "))
		gL = collect(gL[1]:gL[2])
	end
	if cluster == "local"
		addprocs(workers)
	else
		@eval using ClusterManagers
		if cluster == "slurm"
			@eval ClusterManagers.addprocs_slurm(workers)
		elseif cluster == "sge"
			@eval ClusterManagers.addprocs_sge(workers)
		elseif cluster == "htcondor"
			@eval ClusterManagers.addprocs_htc(workers)
		end
	end

	@eval @everywhere using Analytical
	@eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac,rho=$rho,thetaMidNeutral=$theta,al=$shape)

	@eval convolutedSamples = Analytical.binomialDict()
	@eval Analytical.binomOp!($adap,$convolutedSamples.bn);
	@time @eval df = Analytical.rates(param = $adap,convolutedSamples=$convolutedSamples,gH=collect($tmpStrong[1]:$tmpStrong[2]),gL=$gL,gamNeg=collect($tmpNeg[1]:$tmpNeg[2]),iterations = $solutions,shape=$adap.al,output=$output);
end

"Summary statistics from rates. Please provide an analysis folder containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"
@main function summaries(;ne::Int64=1000, samples::Int64=500,rates::String="rates.jld2",summstatSize::Int64=100000,replicas::Int64=100,analysis::String="<folder>",bootstrap::String="true",workers::Int64=1)
	
	#=ne=1000;samples=661;rates="/home/jmurga/ratesBgs.jld2";summstatSize=100000;replicas=100;analysis="/home/jmurga/mkt/202004/rawData/summStat/tgp/wg";bootstrap="true"=#

	addprocs(workers)
	
	@eval @everywhere using Analytical, DataFrames, CSV, PoissonRandom, ProgressMeter
	@eval using JLD2
	
	addprocs(workers)
	
	@eval @everywhere using Analytical, DataFrames, CSV, PoissonRandom, ProgressMeter
	@eval using JLD2
	
    @eval h5file    = jldopen($rates)

    @eval adap      = Analytical.parameters(N=$ne,n=$samples)
    @eval adap.dac  = h5file[string($ne) * "/" * string($samples) * "/dac"]

    @eval if $bootstrap == "true"
    	@eval summstat  = Analytical.summaryStatsFromRates(;param=$adap,rates=$h5file,analysisFolder=$analysis,summstatSize=$summstatSize,replicas=$replicas,bootstrap=true)
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
