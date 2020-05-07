function readData(file)
	df = CSV.read(file)

	resampling = vcat([10,25,50,75],collect(100:100:1000),collect(1000:500:4000))
	tmp = zeros(length(resampling),2)
	for i in 1:length(resampling)
		idx = StatsBase.sample(axes(df, 1), resampling[i]; replace = true, ordered = true)
		tmp[i,:] = sum(convert(Array, df[idx,:]), dims = 1)
	end

	out = vcat(convert(Matrix,unique(df)),tmp,sum(convert(Array, df), dims = 1))
	return out
end

function parseSfs(data,sfsColumns,divColumns)

	df = CSV.read(data,header=false,delim=' ')

    tmp  = split.(df[:,sfsColumns], ",")
    f(x) = Parsers.parse.(Float64,x[2:end-1])
    pn   = round.(reduce(vcat,tmp[:,1] .|> f),digits=4) |> StatsBase.countmap
    ps   = round.(reduce(vcat,tmp[:,2] .|> f),digits=4) |> StatsBase.countmap

    x = zeros(adap.nn)
	y = zeros(adap.nn)
    for i in 1:adap.nn
        try
            x[i] = pn[round.((i/adap.nn),digits=4)]
            y[i] = ps[round.((i/adap.nn),digits=4)]
        catch
            x[i] = 0
            y[i] = 0
        end
    end

	sfs = x .+ y

	P  = sum(sfs)
	D = convert(Matrix,df[:,divColumns]) |> sum
    return [P,sfs,D]
end

function meanQ(x,column=5)
	m = mean(x[:,column])
	q = Statistics.quantile(x[:,column],[0.05,0.95])
	return append!(q,m)
end

function ABCreg(;data::String, prior::String, nparams::Int64, nsummaries::Int64, outputPath::String, outputPrefix::String,tolerance::Float64, regressionMode::String,regPath="/home/jmurga/ABCreg/src/reg")

	reg = `$regPath -p $prior -d $data -P $nparams -S $nsummaries -b $outputPath/$outputPrefix -$regressionMode -t $tolerance`
	run(reg)

	files = filter(x -> occursin(outputPrefix,x), readdir(outputPath))
	# openFiles(f) = GZip.open(DelimitedFiles.readdlm,outputPath*"/"*f)
	openFiles(f) = convert(Matrix,CSV.read(GZip.open(outputPath*"/"*f),header=false))

	estimates = Array{Float64}(undef,length(files),3)
	estimates = files .|> openFiles .|> meanQ
	results = reduce(hcat,estimates) |> transpose
	results = convert(Matrix,results)
	# results = Dict(zip(files, estimates))

	rm.(outputPath.*"/".*files)

	return results
end
