using Fire, Distributed

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

"Summary statistics from rates"
@main function summStat(;ne::Int64=1000, samples::Int64=500,rates::String="rates.jld2",summstatSize::Int64=100000,replicas::Int64=100,analysis::String="<folder>",output::String="<output>",bootstrap::Bool=false,workers::Int64=1)
	

	#=ne=1000;samples=500;rates="/home/jmurga/rates.jld2";summstatSize=100000;replicas=100;analysis="/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.1_0.999";output="<output>";bootstrap=false=#
	#=addprocs(workers)=#

	function openSfsDiv(x,y,dac,h=false,bootstrap=false)

		sfs = CSV.read(x,header=h,DataFrame) |> Array

		if bootstrap
			sfs[:,2:end] = PoissonRandom.pois_rand.(sfs[:,2:end])
		end

		scumu = Analytical.cumulativeSfs(sfs)
		s = sum(scumu[:,2:3],dims=2)[dac]
		divergence = CSV.read(y,DataFrame,header=h) |> Array
		d = [sum(divergence[1:2])]
		α = @. 1 - (divergence[2]/divergence[1] * scumu[:,2]/scumu[:,3])[dac]
		α = round.(α,digits=5)
		return(s,d,α)	
	end

	@eval @everywhere using Analytical
	@eval using JLD2, DataFrames, CSV, PoissonRandom
	
	@eval h5file = jldopen($rates)

    @eval adap     = Analytical.parameters(N=$ne,n=$samples)
    @eval adap.dac = h5file[string($ne) * "/" * string($samples) * "/dac"]

	run(`mkdir -p $output`)
	r = collect(1:replicas)
	if bootstrap
		sFile  = fill(analysis * "/sfs.tsv",replicas)
		dFile  = fill(analysis * "/div.tsv",replicas)
		header = fill(true,replicas)
		b      = fill(true,replicas)
	else
		sFile  = @. analysis * "/sfs" * string(r) * ".tsv"
		dFile  = @. analysis * "/div" * string(r) * ".tsv"
		header = fill(false,replicas)
		b      = fill(false,replicas)
	end

	println("here3")
	@eval tmp = openSfsDiv.($sFile,$dFile,fill($adap.dac,$replicas),$header,$b)

	sfs = [tmp[i][1] for i=1:replicas]
	d   = [tmp[i][2] for i=1:replicas]
	α   = [tmp[i][3] for i=1:replicas]

	@eval summstat = Analytical.summaryStatsFromRates(param=$adap,rates=$h5file,divergence=$d,sfs=$sfs,summstatSize=$summstatSize,replicas=$replicas)

	for i in eachindex(α)
		CSV.write(output * "/alpha_" * string(i) * ".tsv",DataFrame(repeat(α[i]',2)),delim='\t',header=false)
		CSV.write(output * "/" * "summstat_" * string(i) * ".tsv",DataFrame(summstat[i]),delim='\t',header=false)
	end
end

"ABCreg inference"
@main function abcInference(;abcreg::String="<path to ABCreg>",parallel::String="<path to GNU parallel>")
	
	addprocs(workers)
end
