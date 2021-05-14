# Summary statistics
To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md)

Before starting the summary statistics, consider parallelizing the process using Julia Distributed computing. If you are following the tutorial step by step, do not input the following commands.

```julia
using Distributed
addprocs(7)

@everywhere using Analytical, CSV, DataFrames, JLD2, ProgressMeter
```

Before to compute the summary statistic, load your analytical rates and declare a model specifying the samples size and a DAC.

```julia
# Opening rates
h5file   = jldopen("analysis/rates.jld2")
adap = Analytical.parameters(n=661,dac=[2,4,5,10,20,50,200,661,925])
```

Note you can only input DAC already estimated, nonetheles you can perform any subset from the estimated DAC. To check the estimated DAC you can follow the hierarchy of the ```h5file``` variable.

```julia 
# Check hierarchy
h5file

JLDFile /home/jmurga/analysis/rates.jld2 (read-only)
 └─📂 1000
    └─📂 661
       ├─🔢 models
       ├─🔢 neut
       ├─🔢 sel
       ├─🔢 dsdn
       └─🔢 dac
```

```julia
# Checking estimated dac, string pattern inside the HDF5 variable
h5file["1000/661/dac"]

14-element Vector{Int64}:
    1
    2
    4
    5
   10
   20
   50
  100
  200
  400
  500
  661
  925
 1000
```

To standardize the summary statistic estimation, the function ```Analytical.summaryStatsFromRates``` will search and read the SFS and divergence files given a folder. Please be sure that you write the SFS and divergence files (check [Parsing genomic data](input.md)) using the prefix *sfs* and *div* to read the files correctly. Otherwise, the function will not read the files correctly.

We include the argument ```bootstrap``` to perform bootstrap analysis following [polyDFE](https://github.com/paula-tataru/polyDFE) manual. In the following example we boostrap the SFS and divegence file 100 times subseting 10^5 summary statistic for each dataset:

```julia
@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="analysis/",summstatSize=10^5,replicas=100,bootstrap=true);
```

The function will create a summary statistic file and the observed data file (*summaries.txt* and *alphas.txt* respectively). Both files will be used to perform the ABC inference. Each line in *alphas.txt* contains the $\alpha_(x)$ estimations from the bootstrapped SFS.
