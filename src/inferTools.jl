

function parseSfs(;data,output,sfsColumns::Array{Int64,1}=[3,5],divColumns::Array{Int64,1}=[6,7])
	
	g(x) = parse.(Float64,x[2:end-1])
	
	if(data isa String)
		P   = Array{Int64}(undef,1)
		D   = Array{Int64}(undef,1)
		sfs = Array{Float64}(undef, adap.nn -1 ,1)
		newData = Array{Float64}(undef, 1,4)

		df = read(data,header=false,delim=' ')

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

		# Saving summarize data to abc. Ds, Dn, Ps, Pn
		newData = [sum(df[:,divColumns[1]]) sum(df[:,divColumns[2]]) sum(y) sum(x)]

		# Empirical data to analytical estimations
		sfs .= x .+ y
		P = sum(sfs)
		D = convert(Matrix,df[:,divColumns]) |> sum
		
		write(output, DataFrame(newData), delim='\t',writeheader=false)
		return [P,sfs,D]

	else
		P   = Array{Int64}(undef,length(data))
		D   = Array{Int64}(undef,length(data))
		sfs = Array{Int64}(undef,length(data),adap.nn-1)
		newData = Array{Int64}(undef,length(data),4)

		for i in 1:length(data)
			df = read(data[i],header=false,delim=' ')
	
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
			
			# Saving summarize data to abc. Ds, Dn, Ps, Pn
			newData[i,:] = [sum(df[:,divColumns[1]]) sum(df[:,divColumns[2]]) sum(y) sum(x)]
			
			# Empirical data to analytical estimations
			sfs[i,:] = x .+ y
			P[i]  = sum(vcat(sfs[i])...)
			D[i] = convert(Matrix,df[:,divColumns]) |> sum

		end

		write(output, DataFrame(newData),delim='\t',writeheader=false)
		return [P,permutedims(sfs),D]

	end
end

function meanQ(x,columns=[5,6,7])
	x = x[:,columns]
	m = StatsBase.mean(x,dims=1)
	
	qt = Array{Float64}(undef,size(x,2),2)
	for i in 1:size(x,2)
		qt[i,:] = StatsBase.quantile(x[:,i],[0.05,0.95])
	end

	return vcat(m,permutedims(qt))
end

function ABCreg(;data::String, prior::String, nparams::Int64, nsummaries::Int64, outputPath::String, outputPrefix::String,tolerance::Float64, regressionMode::String,regPath="/home/jmurga/ABCreg/src/reg")

	reg = `$regPath -p $prior -d $data -P $nparams -S $nsummaries -b $outputPath/$outputPrefix -$regressionMode -t $tolerance`

	openFiles(f) = convert(Matrix,read(open(f),header=false))

	run(reg)

	files = outputPath .* filter(x -> occursin(outputPrefix,x), readdir(outputPath))

	posteriors = files .|> openFiles
	estimates  = posteriors .|> meanQ
	
	return posteriors,estimates
end

function plotPosterior(data,file,imgSize)

	Plots.gr()
	Plots.theme(:wong2)
	
	p1 = StatsPlots.density(data[:,[5,6,7]],legend = :topright, fill=(0, 0.3),xlabel = "alpha",label = ["alpha strong" "alpha weak" "alpha"],ylabel = "Posterior density", lw = 0.5,fmt = :svg,bottom_margin=10mm,left_margin=10mm,size=imgSize)
	Plots.savefig(file)

end

function readData(file)
	df = read(file)

	resampling = vcat([10,25,50,75],collect(100:100:1000),collect(1000:500:4000))
	tmp = zeros(length(resampling),2)
	for i in 1:length(resampling)
		idx = StatsBase.sample(axes(df, 1), resampling[i]; replace = true, ordered = true)
		tmp[i,:] = sum(convert(Array, df[idx,:]), dims = 1)
	end

	out = vcat(convert(Matrix,unique(df)),tmp,sum(convert(Array, df), dims = 1))
	return out
end
