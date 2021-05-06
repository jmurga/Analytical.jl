# Input data

To estimate summary statistics, you need to provide empirical SFS and divergence files. As explained in section [data](data.md), you can directly parse TGP or DGN data using our module. Nonetheless, you can input any other SFS and divergence file.

You can easily use Julia to download the files using Julia or Bash. You also can download the files manually

```julia
download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt","/home/jmurga/tgpData/tgp.txt")
```

The function ```Analytical.parseSfs``` will parse the raw data, creating two variables of type: ``Matrix{Float64}``` and ```Vector{Int64}``` required to estimate summary statistics.

```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "/home/jmurga/tgpData/tgp.txt")
```

To standardize the estimation next steps, search and read the SFS and divergence file automatically. To do that, you need to create an analysis folder. Please write your SFS and divergence data in the analysis folder using the prefix *sfs* and *div*

```julia
CSV.write("/home/jmurga/tgpData/sfsTgp.tsv",DataFrame(sfs),delim='\t',header=false)
CSV.write("/home/jmurga/tgpData/divTgp.tsv",DataFrame(permutedims(divergence)),delim='\t',header=false)
```

Once you have estimated (or download) the analytical rates and parsed the SFS and divergence information, you can estimate the summary statistics.