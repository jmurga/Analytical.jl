using Fire, Distributed

function openSfsDiv(x,y,dac,replicas,bootstrap="false")

	sfs = Array.(CSV.read.(x,DataFrame))
	divergence = Array.(CSV.read.(y,DataFrame))

	if bootstrap == "true"
		sfs = repeat(sfs,replicas)
		pr(x) = hcat(x[:,1],PoissonRandom.pois_rand.(x[:,2:end]))
		sfs = pr.(sfs)
		divergence = repeat(divergence,replicas)

	end

	scumu = Analytical.cumulativeSfs.(sfs)
	f(x,d=dac) = sum(x[:,2:3],dims=2)[d]
	s = f.(scumu)


	d = [[sum(divergence[i])] for i in eachindex(divergence)]
	al(a,b,c=dac) = @. round(1 - (b[2]/b[1] * a[:,2]/a[:,3])[c],digits=5)
	α = al.(scumu,divergence)
	return(s,d,α)	
end

"Function to estimate rates"
@main function rates(;ne::Int64=1000, samples::Int64=500, gamNeg::String="-1000 -200", gL::String="5 10", gH::String="400 1000",dac::String="2,4,5,10,20,50,200,500,700",solutions::Int64=1000000,output::String="/home/jmurga/rates.jld2",workers::Int64=1)

	tmpNeg    = parse.(Int,split(gamNeg," "))
	tmpWeak   = parse.(Int,split(gL," "))
	tmpStrong = parse.(Int,split(gH," "))
	dac       = parse.(Int,split(dac,","))

	addprocs(workers)

	@eval @everywhere using Analytical
	@eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac)

	@eval convolutedSamples = Analytical.binomialDict()
	@eval Analytical.binomOp!($adap,$convolutedSamples.bn);
	@time @eval df = Analytical.ratesToStats(param = $adap,convolutedSamples=$convolutedSamples,gH=collect($tmpStrong[1]:$tmpStrong[2]),gL=collect($tmpWeak[1]:$tmpWeak[2]),gamNeg=collect($tmpNeg[1]:$tmpNeg[2]),iterations = $solutions,shape=$adap.al,output=$output);
end

"Summary statistics from rates. Please provide an analysis folder containing the SFS and divergence file. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/"

@main function summStat(;ne::Int64=1000, samples::Int64=500,rates::String="rates.jld2",summstatSize::Int64=100000,replicas::Int64=100,analysis::String="<folder>",bootstrap::String="true",workers::Int64=1)
	
	#=ne=1000;samples=661;rates="/home/jmurga/ratesBgs.jld2";summstatSize=100000;replicas=100;analysis="/home/jmurga/mkt/202004/rawData/summStat/tgp/wg";bootstrap="true"=#

	addprocs(workers)
	
	@eval @everywhere using Analytical
	@eval using JLD2, DataFrames, CSV, PoissonRandom, ProgressMeter
	
	@eval h5file = jldopen($rates)

	@eval adap     = Analytical.parameters(N=$ne,n=$samples)
	@eval adap.dac = h5file[string($ne) * "/" * string($samples) * "/dac"]

	@eval sFile = filter(x -> occursin("sfs",x), readdir($analysis,join=true))
	@eval dFile = filter(x -> occursin("div",x), readdir($analysis,join=true))

	@eval sfs,d,α = openSfsDiv($sFile,$dFile,$adap.dac,$replicas,$bootstrap)

	@eval summstat = Analytical.summaryStatsFromRates(param=$adap,rates=$h5file,divergence=$d,sfs=$sfs,summstatSize=$summstatSize,replicas=$replicas)

	@eval w(x,name) = CSV.write(name,DataFrame(x),delim='\t',header=false)

	@eval w.(permutedims.($α),@. $analysis * "/alpha_" * string(1:$replicas) * ".tsv")
	@eval w.(summstat,@. $analysis * "/summstat_" * string(1:$replicas) * ".tsv")
	#=
	for i in eachindex(α)
		@eval CSV.write($analysis * "/alpha_" * string($i) * ".tsv",DataFrame($α[$i]'),delim='\t',header=false)
		@eval CSV.write($analysis * "/" * "summstat_" * string(i) * ".tsv",DataFrame($summstat[$i]),delim='\t',header=false)
	end=#
end

"ABCreg inference"
@main function abcInference(;analysis::String="<folder>",replicas::Int64=100,P::Int64=5,S::Int64=9,tol::Float64=0.01,workers::Int64=1,parallel::String="true")
	
	reg = chomp(read(`which reg`,String))

	@eval aFile = @. $analysis * "/alpha_" * string(1:$replicas) * ".tsv"
	@eval sumFile = @. $analysis * "/summstat_" * string(1:$replicas) * ".tsv"
	@eval out = @. $analysis * "/out_" * string(1:$replicas)
	if parallel == "true"
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

		out = filter(x -> occursin("post",x), readdir(analysis,join=true))
		out = filter(x -> !occursin(".1.",x),out)

		open(x) = Array(CSV.read(GZip.open(x),DataFrame,header=false))
		@eval posteriors = open.($out)
		
		@eval getmap(x) = rcopy(R"""matrix(apply(x,2,getmap),nrow=1)""")
		@eval tmp = getmap.($posteriors)
		@eval maxp = DataFrame($tmp,[:aw,:as,:a,:gamNeg,:shape])
	
		al  = maxp[:,1:3]
		gam  = maxp[:,4:end]

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
