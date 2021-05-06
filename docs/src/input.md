# Input data

To estimate summary statistics you need to provide empirical SFS and divergence files. As explained in section [Data](@ref), you can directly parse TGP or DGN data using our module. Nonetheless you ca input any other SFS and divergence file.

You can easly use julia to download the files using Julia or Bash. You also can download the files manually


```julia
download("https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt","/home/jmurga/tgpData/tgp.txt")
```

Parsing the data to SFS and divergence. The function ```Analytical.parse``` will create two variable of type ```Matrix{Float64}``` and ```Vector{Int64}``` required to estimate summary statistics.

```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "/home/jmurga/tgpData/tgp.txt")
```

To standarize the estimation next steps search and read automatically the SFS and divergence file. To do that, you need to create an analysis folder. Please write your SFS and divergence data in the analysis folder using the prefix *sfs* and *div*

```julia
CSV.write("/home/jmurga/tgpData/sfsTgp.tsv",DataFrame(sfs),delim='\t',header=false)
CSV.write("/home/jmurga/tgpData/divTgp.tsv",DataFrame(permutedims(divergence)),delim='\t',header=false)
```


Once you have estimated (or download) the analytical rates and you have parsed the SFS and divergence information, you can estimate the summary statistics.