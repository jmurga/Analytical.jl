# Summary statistics
To estimate summary statistics we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](@ref)

Before to start the summary statistics, consider parallelize the process using Julia Distributed computing. If you are following the tutorial step by step, do not input the next commands

```julia
using Distributed
addprocs(7)

@everywhere using Analytical, CSV, DataFrames, JLD2, ProgressMeter
```

Load your analytical rates and declare a model specifing a dac. The dac will be used to subset summary statistics. You only can input dac values already estimated. To check the selected dac at rates estimations, you can access to the H5 file

```julia
# Opening rates
h5file   = jldopen("/home/jmurga/rates.jld2")
# Checking estimated dac
adap = Analytical.parameters(n=661)
h5file["1000/" * string(adap.n) * "/dac"]
# Selecting dac to perform summary statistics
adap.dac = [2,4,5,10,20,50,200,661,925]
```

To standarize the summary statistic estimation, the function ```Analytical.summaryStatsFromRates``` will search and read the SFS and divergence files into an analysis folder. To properly read the files, please be sure that you write the SFS and divergence files using the prefix *sfs* and *div*, otherwise, the function will not read the files properly

We include the argument bootstrap to perform boostrap analysis following [polyDFE](https://github.com/paula-tataru/polyDFE)

```julia
@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="/home/jmurga/tgpData/",summstatSize=10^5,replicas=100,bootstrap=true);
```

The function will create summary statistics files and empirical data input at ABC inference