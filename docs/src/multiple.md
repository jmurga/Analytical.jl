# Infering the adapation rate from multiple datasets
As we explained at section [Input Data](input.data) it is possible to subset the TGP dataset given a matrix of Ensembl/Flybase ids. In this section we will use whole-genome protein-coding information, Viral Interacting Proteins and Non-Viral Interacting Proteins to infer the empirical adaptation.

The next step assume that you execute the previous sections. Please before to execute this example, be sure you already performed the [rates estimates](rates.md) and your *analysis/* folder contain the TGP data.

We are going to create three separate directories into the *analysis/* directory we are working with to perform the estimations. You can create the folder manually, through bash, or directly in the Julia interpreter

```julia
using Distributed
addprocs(7)
@everywhere using Analytical
using CSV, DataFrames, JLD2

run(`mkdir -p analysis/wg/ analysis/vips/ analysis/nonvips/`)
```

Once you have the folder please download from our repository the files containing the VIPs and Non-VIPs ids.

```julia
run(`curl -o analysis/vips/vipList.txt https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/vipList.txt`)
run(`curl -o analysis/nonvips/nonvipList.txt https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/nonvipList.txt`)
```

Now we are going to parse our three dataset to output the SFS and divergence files into each analysis folder. We we follow the example provided at [Input Data](input.data) section. 

 - Whole-Genome dataset
```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt")
CSV.write("analysis/wg/sfsTgp.tsv",DataFrame(sfs,:auto),delim='\t',header=false)
CSV.write("analysis/wg/divTgp.tsv",DataFrame(permutedims(divergence),:auto),delim='\t',header=false)
```

 - VIPs dataset
```julia
vipsList = CSV.read("analysis/vips/vipList.txt",DataFrame) |> Array

alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt",geneList=vipsList)
CSV.write("analysis/vips/sfsTgp.tsv",DataFrame(sfs,:auto),delim='\t',header=false)
CSV.write("analysis/vips/divTgp.tsv",DataFrame(permutedims(divergence),:auto),delim='\t',header=false)
```

 - Non-VIPs dataset
```julia
nonvipsList = CSV.read("analysis/nonvips/nonvipList.txt",DataFrame) |> Array

alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt",geneList=nonvipsList)
CSV.write("analysis/nonvips/sfsTgp.tsv",DataFrame(sfs,:auto),delim='\t',header=false)
CSV.write("analysis/nonvips/divTgp.tsv",DataFrame(permutedims(divergence),:auto),delim='\t',header=false)
```

Once you have the data parsed, you can estimate the summary statistic following [Summary statistic](summstat.md) section. Load the rates and select a DAC before to continue the analysis

```julia
h5file   = jldopen("analysis/rates.jld2")
adap = Analytical.parameters(n=661,dac=[2,4,5,10,20,50,200,661,925])
```
 - Whole-Genome dataset
```julia
@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="analysis/wg/",summstatSize=10^5,replicas=100,bootstrap=true);
```

 - VIPs dataset
```julia
@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="analysis/vips/",summstatSize=10^5,replicas=100,bootstrap=true);
```

 - Non-VIPs dataset
```julia
@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="analysis/nonvips/",summstatSize=10^5,replicas=100,bootstrap=true);
```
