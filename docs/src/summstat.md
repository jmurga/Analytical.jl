# Summary statistics
To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md)

Before starting the summary statistics, consider parallelizing the process using Julia Distributed computing. If you are following the tutorial step by step, do not input the following commands.

```julia
using Distributed
addprocs(7)

@everywhere using Analytical, CSV, DataFrames, JLD2, ProgressMeter
```

Load your analytical rates and declare a model specifying a DAC. The selected DAC will be used to subset summary statistics. You only can input DAC values already estimated. To check the selected DAC at rates estimations, you can access the HDF5 file

```julia
# Opening rates
h5file   = jldopen("/home/jmurga/rates.jld2")
# Checking estimated dac
adap = Analytical.parameters(n=661)
h5file["1000/" * string(adap.n) * "/dac"]
# Selecting dac to perform summary statistics
adap.dac = [2,4,5,10,20,50,200,661,925]
```

To standardize the summary statistic estimation, the function ```Analytical.summaryStatsFromRates``` will search and read the SFS and divergence files into an analysis folder. Please be sure that you write the SFS and divergence files using the prefix *sfs* and *div* to read the files correctly. Otherwise, the function will not read the files correctly.

We include the argument ```bootstrap``` to perform bootstrap analysis following [polyDFE](https://github.com/paula-tataru/polyDFE)

```julia
@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="/home/jmurga/tgpData/",summstatSize=10^5,replicas=100,bootstrap=true);
```

The function will create summary statistics files and empirical data input at ABC inference