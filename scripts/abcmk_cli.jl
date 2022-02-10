using ArgParse, Distributed


function parse_cli(args)
    settings = ArgParseSettings("ABC-MK CLI")

    @add_arg_table! settings begin
        # CMD 1
        "rates", "R"
            help = "Function to solve fixation and polymorphic rates analitically. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.\nIf rho and/or theta are set to nothing, the function will input random values given the range 0.0005:0.0005:0.01. Otherwise you can fix the values.\nIf gL is set to nothing, the function will not account the role of the weakly selected alleles in the estimation.\nThe function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.\nPlease check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/Analytical.jl/dev/cli/ to check model
            "
            action = :command
        # CMD 2
        "parse_data", "P"
            help = "Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics. Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli/"
            action = :command
        # CMD 3
        "summaries", "S"
            help = "Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics. Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli/"
            action = :command
        # CMD 4
        "inference", "I"
            help = "ABCreg inference. The function returns posterior distributions from ABC inference. Each posterior file contains information about alpha_w, alpha_s, alpha, gamNeg and shape parameter. The number of posterior distributions will depend on the number of bootstrap replicas. Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli"
            action = :command
        
    end

    add_arg_table!(settings["rates"],
        ["--pop_size"],
        Dict(
            :help => "Population size",
            :arg_type => Int64,
            :default => 1000
        ),
        ["--sample_size"],
        Dict(
            :help => "Sample size",
            :arg_type => Int64,
            :required => true
        ),
        ["--dac"],
        Dict(
            :help => "Derived Allele Count",
            :arg_type => String,
            :required => true
        ),
        ["--gam_neg"],
        Dict(
            :help => "Selection coefficient for deleterious alleles",
            :arg_type => String,
            :required => true
        ),
        ["--positive_strong"],
        Dict(
            :help => "Selection coefficient for strongly beneficial alleles",
            :arg_type => String,
            :required => true
        ),
        ["--positive_weak"],
        Dict(
            :help => "Selection coefficient for weakly beneficial alleles",
            :arg_type => String,
            :required => true
        ),
        ["--shape"],
        Dict(
            :help => "Shape value modeling Gamma distribution for deleterious alleles",
            :arg_type => Float64,
            :default => 0.184
        ),
        ["--rho"],
        Dict(
            :help => "Recombination rate",
            :arg_type => Float64,
            :default => 0.001
        ),
        ["--theta"],
        Dict(
            :help => "Mutation rate on coding locus",
            :arg_type => Float64,
            :default => 0.001
        ),
        ["--solutions"],
        Dict(
            :help => "Mutation rate on coding locus",
            :arg_type => Int64,
            :default => 100_000
        ),
        ["--output"],
        Dict(
            :help => "Output file",
            :arg_type => String,
            :default => "rates.jld2"
        ),
        ["--scheduler"],
        Dict(
            :help => "Select scheduler manager",
            :arg_type => String,
            :default => "local"
        ),
        ["--nthreads"],
        Dict(
            :help => "Select number of threads to parallelize",
            :arg_type => Int64,
            :default => 1
        ),
    )

    add_arg_table!(settings["parse_data"],
        "folder",
        Dict(
            :help => "a positional argument",
            :required => true,
            :arg_type => String
        ),
        ["--dataset","-d"],
        Dict(
            :help => "a positional argument",
            :default => "tgp",
            :arg_type => String
        ),
        ["--gene_list","-g"],
        Dict(
            :help => "a positional argument",
            :arg_type => Union{Bool,String},
            :default => false
        ),
        ["--bins","-b"],
        Dict(
            :help => "a positional argument",
            :arg_type => Union{Bool,Int64},
            :default => false
        )
    )

    add_arg_table!(settings["summaries"],
        "folder",
        Dict(
            :help => "Folder path containing SFS and divergence files to run the analysis",
            :required => true,
            :arg_type => String
        ),
        ["--rates"],
        Dict(
            :help => "H5 file containing precomputed rates",
            :required => true,
            :arg_type => String
        ),       
        ["--sample_size"],
        Dict(
            :help => "Sample size",
            :required => true,
            :arg_type => Int64
        ),
        ["--dac"],
        Dict(
            :help => "Derived allele count",
            :required => true,
            :arg_type => String
        ),
        ["--summstat_size"],
        Dict(
            :help => "Define number of summary estatistic to perform ABC",
            :required => true,
            :arg_type => Int64
        ),
        ["--bootstrap"],
        Dict(
            :help => "Allow bootstrap following polyDFE manual",
            :arg_type => Bool,
            :default => false
        ),
        ["--replicas"],
        Dict(
            :help => "Number of bootstrap replicas",
            :arg_type => Int64,
            :default => 1
        ),
        ["--nthreads"],
        Dict(
            :help => "Select scheduler manager",
            :arg_type => Int64,
            :default => 1
        ),
        ["--scheduler"],
        Dict(
            :help => "Select scheduler manager",
            :arg_type => String,
            :default => "local"
        )
    )
    add_arg_table!(settings["inference"],
        "folder",
        Dict(
            :help => "Folder path containing SFS and divergence files to run the analysis",
            :required => true,
            :arg_type => String
        ),
        ["--S"],
        Dict(
            :help => "Define number of summary estatistic to perform ABC",
            :required => true,
            :arg_type => Int64
        ),
        ["--tol"],
        Dict(
            :help => "Tolerance",
            :required => true,
            :arg_type => Float64
        ),
        ["--abcreg"],
        Dict(
            :help => "ABCreg static binary",
            :required => true,
            :arg_type => String
        ),
        ["--nthreads"],
        Dict(
            :help => "Select scheduler manager",
            :arg_type => Int64,
            :default => 1
        ),
        ["--scheduler"],
        Dict(
            :help => "Select scheduler manager",
            :arg_type => String,
            :default => "local"
        )
    )

    return parse_args(settings)
end

cli = parse_cli(ARGS)

for cmd in keys(cli)
    if(cmd == "parse_data")

        @eval using Analytical, DataFrames, CSV
        
        folder = cli[cmd]["folder"]

        run(`mkdir -p $folder`)

        dataset = lowercase(cli[cmd]["dataset"])
        data    = folder * "/" * dataset * ".txt"

        @eval download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/"* $dataset * ".txt",$data)

        # Check if bins or gene_list are defined
        gene_list = cli[cmd]["gene_list"]
        @eval if $gene_list != false
            @eval gList = CSV.read($gene_list,DataFrame,header=false) |> Array
        else
            @eval gList = nothing
        end

        bins = cli[cmd]["bins"]

        @eval if $bins != 0
            @eval bins_size = $bins
        else
            bins_size = nothing
        end

        # Parsing TGP data
        if dataset == "tgp"
            @eval α,sfs, divergence = Analytical.parse_sfs(sample_size=661,data=$data,gene_list=$gList,bins=$bins_size)
        elseif occursin("zi",dataset)
            @eval α,sfs, divergence = Analytical.parse_sfs(sample_size=154,data=$data,gene_list=$gList,bins=$bins_size,isolines=true)
        elseif occursin("ral",dataset)
            @eval α,sfs, divergence = Analytical.parse_sfs(sample_size=160,data=$data,gene_list=$gList,bins=$bins_size,isolines=true)
        end
        # Writting data to folder
        @eval sName = $folder * "/sfs.tsv"
        @eval dName = $folder * "/div.tsv"

        @eval CSV.write($sName,DataFrame($sfs,:auto),delim='\t',header=false)
        @eval CSV.write($dName,DataFrame($divergence',:auto),delim='\t',header=false)
    elseif (cmd == "rates")

        neg    = parse.(Int,split(cli[cmd]["gam_neg"],":"))
        strong = parse.(Int,split(cli[cmd]["positive_strong"],":"))
        
        if (cli[cmd]["positive_weak"] == false)
            weak = nothing
        else
            weak = parse.(Int,split(cli[cmd]["positive_weak"],":"))
        end

        dac = parse.(Int,split(cli[cmd]["dac"],","))


        if (cli[cmd]["rho"] == false)
            rho = nothing
        else
            rho = cli[cmd]["rho"]
        end

        if (cli[cmd]["theta"] == "nothing")
            theta = nothing
        else
            theta = cli[cmd]["theta"]
        end

        scheduler = cli[cmd]["scheduler"];nthreads = cli[cmd]["nthreads"]
        if scheduler == "slurm"
            @eval using ClusterManagers
            @eval addprocs_slurm($nthreads)
        elseif scheduler == "htcondor"
            @eval using ClusterManagers
            @eval addprocs_htc($nthreads)
        else
            @eval addprocs($nthreads)
        end
        
        ne = cli[cmd]["pop_size"];samples= cli[cmd]["sample_size"];shape = cli[cmd]["shape"]
        @eval @everywhere using Analytical, ParallelUtilities
        @eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac,al=$shape)

        @eval cnv = Analytical.binomial_dict()
        @eval Analytical.binomOp!($adap,$cnv.bn);

        solutions = cli[cmd]["solutions"]; output = cli[cmd]["output"]
        @eval Analytical.rates(param = $adap,convoluted_samples=$cnv,gH=$strong[1]:$strong[2],gL=$weak[1]:$weak[2],gamNeg=$neg[1]:$neg[2],iterations = $solutions,rho=$rho,theta=$theta,shape=$adap.al,output=$output,scheduler=$scheduler);

        for i in workers()
            rmprocs(i)
        end
    elseif (cmd == "summaries")

        scheduler = cli[cmd]["scheduler"];nthreads = cli[cmd]["nthreads"];

        if scheduler == "slurm"
            @eval using ClusterManagers
            @eval addprocs_slurm($nthreads)
        elseif scheduler == "htcondor"
            @eval using ClusterManagers
            @eval addprocs_htc($nthreads)
        else
            @eval addprocs($nthreads)
        end

        samples       = cli[cmd]["sample_size"]
        dac           = parse.(Int,split(cli[cmd]["dac"],","))
        bootstrap     = cli[cmd]["bootstrap"]
        replicas      = cli[cmd]["replicas"]
        rates         = cli[cmd]["rates"]
        summstat_size = cli[cmd]["summstat_size"]
        folder        = cli[cmd]["folder"]

        @eval  using JLD2, DataFrames, CSV, ProgressMeter
        @eval @everywhere using  Analytical, ParallelUtilities
        @eval adap    = Analytical.parameters(n=$samples,dac = $dac)

        @eval if ($bootstrap == true)
            @eval summstat  = Analytical.summary_statistics(param=$adap,h5_file=$rates,analysis_folder=$folder,summstat_size=$summstat_size,replicas=$replicas,bootstrap=true)
        else
            @eval summstat  = Analytical.summary_statistics(param=$adap,h5_file=$rates,analysis_folder=$folder,summstat_size=$summstat_size,replicas=$replicas,bootstrap=false)
        end

        for i in workers()
            rmprocs(i)
        end
    elseif (cmd == "inference")

        scheduler = cli[cmd]["scheduler"]
        nthreads = cli[cmd]["nthreads"]
        folder   = cli[cmd]["folder"]
        S        = cli[cmd]["S"]
        tol      = cli[cmd]["tol"]
        abcreg   = cli[cmd]["abcreg"]

        if scheduler == "slurm"
            @eval using ClusterManagers
            @eval addprocs_slurm($nthreads)
        elseif scheduler == "htcondor"
            @eval using ClusterManagers
            @eval addprocs_htc($nthreads)
        else
            @eval addprocs($nthreads)
        end
        
        @eval @everywhere using Analytical,ParallelUtilities

        @eval Analytical.ABCreg(analysis_folder=$folder,S=$S,tol=$tol,abcreg=$abcreg)

        for i in workers()
            rmprocs(i)
        end
    end
end
