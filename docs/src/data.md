# Parsing genomic data
The module includes functions to parse TGP from Uricchio et al. (2019) and DGN from Murga-Moreno et al. (2019). In addition, the module have a function to parse SFS and divergence from multi-FASTA data following Murga-Moreno et al. (2019)

Please to parse raw data into SFS and divergence counts, first download raw files deposited in our repository:  

 - [TGP](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt)
 - [DGN Zambia population](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/dgnRal.txt)  
 - [DGN Raleigh population](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/dgnZi.txt)  

```bash
mkdir -p analysis/
curl -o analysis/tgp.txt https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt
curl -o analysis/dgnRal.txt https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/dgnRal.txt
```
## Parsing TGP and DGN data manually
Once you have downloaded the files, you can use the function ```Analytical.parseSfs``` to convert the data into SFS and divergence counts. Please check [`Analytical.parseSfs`](@ref) to get more info o execute:

```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt")
```

To save the data, you can use CSV and DataFrames packages

```julia
using CSV, DataFrames
CSV.write("analysis/tgpSfs.tsv",DataFrame(sfs,:auto),delim='\t',header=false)
CSV.write("analysis/tgpDiv.tsv",DataFrame(permutedims(divergence),:auto),delim='\t',header=false)
```

It is possible to directly subset genes IDs using Ensembl or Flybase id. Use a variable of type ```Matrix{String}``` into the argument *geneList*

```julia
download('https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/ensemblList.txt','analysis/ensemblList.txt')
ensemblList = CSV.read("analysis/ensemblList.txt",header=false,DataFrame) |> Array

alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt",geneList = ensemblList)
```

If you are going to parse DGN, you need to change the value of the argument *isoline* to *true*. Following the Murga-Moreno et al. (2019) sample size for each population is:

 - Zambia population: 154
 - RAL population: 160

```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 160, data = "analysis/dgnRal.txt",isolines=true)
```