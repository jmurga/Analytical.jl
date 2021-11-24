# Input data

To estimate summary statistics, you need to provide empirical SFS and divergence files. As explained in section [data](data.md), you can directly parse TGP or DGN data using our module. Nonetheless, you can input any other SFS and divergence file.

You can easily use Julia to download the files using Julia or Bash. You also can download the files manually.

```julia
download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt","analysis/tgp.txt")
```

To standardize the estimation, the summary statitics ABC inference functions will search SFS, divergence files into a folder containing both file with prefix *sfs* and *div*. The function ```Analytical.parseSfs``` will parse the raw data, creating two variables of type: ```Matrix{Float64}``` and ```Vector{Int64}``` required to estimate summary statistics.

```julia
alpha, sfs, divergence = Analytical.parse_sfs(sample_size = 661, data = "analysis/tgp.txt")
```

```julia
using CSV, DataFrames

CSV.write("analysis/sfs_tgp.tsv",DataFrame(sfs,:auto),delim='\t',header=false)
CSV.write("analysis/div_tgp.tsv",DataFrame(permutedims(divergence),:auto),delim='\t',header=false)
```

Once you have estimated (or download) the analytical rates and parsed the SFS and divergence information, you can estimate the summary statistics.