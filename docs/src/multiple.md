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

Now you can perform the ABC inference using ABCreg or another ABC software using the files *alphas.txt* and *summaries.txt* deposited in each folder. We will perform the inference using ABCreg. Please compile ABCreg before to perform the execution as described in the [Installation](index.md)

 - Whole-Genome dataset
```julia
Analytical.ABCreg(analysisFolder="analysis/wg/",P=5,S=size(adap.dac,1),tol=0.0001,abcreg="/home/jmurga/ABCreg/src/reg");
```

 - VIPs dataset
```julia
Analytical.ABCreg(analysisFolder="analysis/vips/",P=5,S=size(adap.dac,1),tol=0.0001,abcreg="/home/jmurga/ABCreg/src/reg");
```

 - Non-VIPs dataset
```julia
Analytical.ABCreg(analysisFolder="analysis/nonvips/",P=5,S=size(adap.dac,1),tol=0.0001,abcreg="/home/jmurga/ABCreg/src/reg");
```

Once the posterior distributions are estimated you can perform the Maximum-A-Posterior (MAP) estimates. We performed the MAP estimates following ABCreg examples. Please be sure you have installed R and the packages(ggplot2, locfit and data.table). This packages are already installed in Docker and Singularity images. Load the R functions to estimate and plot MAP executing the following command

```julia
Analytical.sourcePlotMapR(script="analysis/script.jl")
```

 - Whole-Genome dataset
```julia
tgpmap = Analytical.plotMap(analysisFolder="analysis/wg/")
```

 - VIPs dataset
```julia
vipsmap = Analytical.plotMap(analysisFolder="analysis/vips/")
```

 - Non-VIPs dataset
```julia
nonvipsmap = Analytical.plotMap(analysisFolder="analysis/nonvips/")
```
