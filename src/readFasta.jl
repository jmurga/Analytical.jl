function sequencesToMatrix(samples::Int64,length::Int64,sequences::Array{Tuple{String,String},1})

	matrix = Array{Char}(undef,samples + 2,length)

	deleteIndex = Array{Int8}[]

	for i in 1:samples
		# Extract each sample sequence
		# tmp = multiFasta[samples[i]][:].seq
		tmp = sequences[i][2]
		matrix[i+1,:] = collect(tmp)
		# if(len(tmp) != seqLen)
		# 	print('errorAlign')
		# 	break
		# if('N' in tmp)
		# 	deleteIndex.append(i)
		# else
	end
	# matrix = np.delete(matrix,deleteIndex,0)
	return matrix
end

function degeneracy(data::String,codonDict::String)
	#degeneracy DICTIONARIES
	standardDict = Dict{String,String}(
		"TTT"=> "002", "TTC"=> "002", "TTA"=> "202", "TTG"=> "202",
		"TCT"=> "004", "TCC"=> "004", "TCA"=> "004", "TCG"=> "004",
		"TAT"=> "002", "TAC"=> "002", "TAA"=> "022", "TAG"=> "002",
		"TGT"=> "002", "TGC"=> "002", "TGA"=> "020", "TGG"=> "000",
		"CTT"=> "004", "CTC"=> "004", "CTA"=> "204", "CTG"=> "204",
		"CCT"=> "004", "CCC"=> "004", "CCA"=> "004", "CCG"=> "004",
		"CAT"=> "002", "CAC"=> "002", "CAA"=> "002", "CAG"=> "002",
		"CGT"=> "004", "CGC"=> "004", "CGA"=> "204", "CGG"=> "204",
		"ATT"=> "003", "ATC"=> "003", "ATA"=> "003", "ATG"=> "000",
		"ACT"=> "004", "ACC"=> "004", "ACA"=> "004", "ACG"=> "004",
		"AAT"=> "002", "AAC"=> "002", "AAA"=> "002", "AAG"=> "002",
		"AGT"=> "002", "AGC"=> "002", "AGA"=> "202", "AGG"=> "202",
		"GTT"=> "004", "GTC"=> "004", "GTA"=> "004", "GTG"=> "004",
		"GCT"=> "004", "GCC"=> "004", "GCA"=> "004", "GCG"=> "004",
		"GAT"=> "002", "GAC"=> "002", "GAA"=> "002", "GAG"=> "002",
		"GGT"=> "004", "GGC"=> "004", "GGA"=> "004", "GGG"=> "004")

	if (codonDict == "standard")
		degenerateCodonTable = standardDict
	end

	degen="";
	vectorIter = collect(1:3:length(data));
	@inbounds for i in vectorIter
		codon = data[i:i+2]

		if ('N' in codon) || ('-' in codon)
			degen = string(degen,codon)
		else
		degen = degen * degenerateCodonTable[codon]
		end
	end

	return degen
end

function iterFastaMatrix(sequenceMatrix::Array{Char,2})

	# output = DataFrame([Float64,Int64,String],[:daf,:div,:degen])
	output = Array{Any,2}[]

	for n in eachcol(sequenceMatrix)

		degen   = n[1]
		aa      = n[end]
		pol     = n[2:end-1]

		# Undefined Ancestra Allele. Try to clean out of the loop
		if (aa == 'N' || aa == '-')
			continue
		# Monomorphic sites. Try to clean out of the loop
		elseif (size(unique(n[2:end][n[2:end] .!= 'N']),1) == 1)
			continue
		else

			if degen == '4'
				functionalClass = "4fold"
			else
				functionalClass = "0fold"
			end


			# Check if pol != AA and monomorphic
			if (size(unique(pol),1) == 1 && unique(pol)!= aa)
					div = 1; af = 0
					tmp = [af div functionalClass]
					push!(output,tmp);
			else
					div    = 0
					an     = size(pol,1)
					ac     = countmap(pol)

					if (aa ∉ keys(ac))
						continue
					else
						da = collect(keys(ac))
						da = da[da .!= aa ]

						if (size(da,1) < 2)
							af = ac[da[1]]/an
						else
							continue
						end
					end
					tmp = [af div functionalClass]
					# println(tmp);
					push!(output,tmp);
			end
		end
	end

	return reduce(vcat,output)
end

function formatSfs(rawSfsOutput::Array{Any,2},samples::Int64,bins::Int64)

	freq = OrderedDict(round.(collect(1:samples)/samples,digits=5) .=> 0)

	sfs = rawSfsOutput[rawSfsOutput[:,2].!=1,[1,3]]
	pn = sort!(OrderedDict(StatsBase.countmap(sfs[sfs[:,2] .!= "4fold",1])))
	ps = sort!(OrderedDict(StatsBase.countmap(sfs[sfs[:,2] .== "4fold",1])))

	sfsPn        = reduceSfs(reduce(vcat,values(merge(+,freq,pn))),bins)'[:,1]
	sfsPs        = reduceSfs(reduce(vcat,values(merge(+,freq,ps))),bins)'[:,1]
	daf          = DataFrame(f=collect(1:bins)/bins,pi=sfsPn,p0=sfsPs)

	divergence = rawSfsOutput[rawSfsOutput[:,2].==1,2:3]

	div = DataFrame(countmap(divergence[:,2]))
	return(daf,div)
end

function uSfsFromFasta(;file::String,reference::String,outgroup::String,samples::Int64,bins::Int64,codonTable::String)

	multiFasta = readfasta(file);
	ref        = readfasta(reference)
	out        = readfasta(outgroup)

	s = samples
	seqLen  = length(ref[1][2])

	multiFastaMatrix = sequencesToMatrix(samples,seqLen,multiFasta);

	degenCode = collect(degeneracy(ref[1][2],codonTable));
	multiFastaMatrix[1,:] = degenCode
	multiFastaMatrix[end,:] = collect(out[1][2])

	m = multiFastaMatrix[:,(multiFastaMatrix[end,:].!= 'N') .& (multiFastaMatrix[end,:].!= '-')]
	m = m[:,(m[1,:] .== '4') .| (m[1,:] .=='0')]

	joinedSfsDiv = iterFastaMatrix(m)

	sfs, div = formatSfs(joinedSfsDiv,samples,bins)
	return sfs,div
end

function uSfsFromVcf(file::String)

	reader = VCF.Reader(GzipDecompressorStream(open(file, "r")))
	output = []

	for record in reader
		ref =  VCF.ref(record)
		alt = VCF.alt(record)
		info = Dict(VCF.info(record))
		aa = info["AA"]
		ac = parse(Int64,info["AC"]);
		an = parse(Int64,info["AN"]);
		af = ac/an

		if (aa == ref)
			daf = af
			tmp = [parse(Int64,VCF.chrom(record)),VCF.pos(record),daf]
			push!(output,tmp)
		elseif (aa == alt)
			daf = 1 - af
			tmp = [parse(Int64,VCF.chrom(record)),VCF.pos(record),daf]
			push!(output,tmp)
		else
			continue
		end
	end
	close(reader)
	return output
end
