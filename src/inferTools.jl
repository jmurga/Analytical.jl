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

function parseSfs(data,sfsColumns=[3,5],divColumns=[5,6])

	g(x) = Parsers.parse.(Float64,x[2:end-1])
	
	if(data isa String)
		P   = Array{Int64}(undef,1)
		D   = Array{Int64}(undef,1)
		sfs = Array{Array}(undef,1)	

		df = CSV.read(data,header=false,delim=' ')

		tmp  = split.(df[:,sfsColumns], ",")
		pn   = round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> StatsBase.countmap
		ps   = round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> StatsBase.countmap

		x = zeros(adap.nn -1)
		y = zeros(adap.nn -1)
		for i in 1:adap.nn -1
			try
				x[i] = pn[round.((i/adap.nn),digits=4)]
				y[i] = ps[round.((i/adap.nn),digits=4)]
			catch
				x[i] = 0
				y[i] = 0
			end
		end

		sfs[i] = x .+ y

		P[i]  = sum(vcat(sfs[i])...)
		D[i]= convert(Matrix,df[:,divColumns]) |> sum

	else
		P   = Array{Int64}(undef,length(data))
		D   = Array{Int64}(undef,length(data))
		sfs = Array{Int64}(undef,length(data),adap.nn-1)
		for i in 1:length(data)
			df = CSV.read(data[i],header=false,delim=' ')
	
			tmp  = split.(df[:,sfsColumns], ",")
			pn   = round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> StatsBase.countmap
			ps   = round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> StatsBase.countmap
	
			x = zeros(adap.nn -1)
			y = zeros(adap.nn -1)
			for i in 1:adap.nn -1
				try
					x[i] = pn[round.((i/adap.nn),digits=4)]
					y[i] = ps[round.((i/adap.nn),digits=4)]
				catch
					x[i] = 0
					y[i] = 0
				end
			end
	
			sfs[i,:] = x .+ y
	
			P[i]  = sum(vcat(sfs[i])...)
			D[i]= convert(Matrix,df[:,divColumns]) |> sum
		end
	end
	return [P,permutedims(sfs),D]
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
