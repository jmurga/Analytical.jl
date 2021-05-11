# Estimating fixation and polymorphic rates considering generalized model of selection and linkage

Before executing the rate estimation, you need to load *Distributed* module and add some threads

```julia
using Distributed
addprocs(7)
```

Then you need to declare the ```Analytical``` module in all the threads using ```@everywhere``` macro. Otherwise, the ```Analytical``` module will perform the estimation just using the main core

```julia
@everywhere using Analytical, ParallelUtilites
using CSV, DataFrames, JLD2, ProgressMeter
```

Declare a variable containing some basic information about your model. We used a sample size of 661 to perform later analysis over TGP data. The selected Derived Alleles Counts (dac) will be used to compute summary statistics and perform ABC inference

```julia
adap              = Analytical.parameters(n=661,dac=[1,2,4,5,10,20,50,100,200,400,500,661,925,1000], al=0.184)
convolutedSamples = Analytical.binomialDict()
Analytical.binomOp!(adap,convolutedSamples.bn)
```

The function ```Analytical.rates``` will perform the analytical estimation of *N* independent models regarding DFE, BGS, mutation rate, and recombination rate. In the following example, we declared the prior distribution for each model parameter. We show how to compute summary statistics in the following section

Note that [```Analytical.rates```](@ref) is the most resource and time-consuming function. In our case, the function will estimate 10^5 independent models. Each model solves the estimation for all possible BGS values. We used BGS values from 0.1 to 0.999 in 5% increments. In total, the example will produce 3.7 million estimates. We have used a hierarchical data structure (HDF5) to facilitate model parameters and rates storage.

The following example took about 1.5 hours to execute on the hardware described at section [Infering the rate and strength of adaptation](empirical.md)

```julia
@time df = Analytical.rates(param = adap,convolutedSamples=convolutedSamples,gH=collect(200:2000),gL=collect(1:10),gamNeg=collect(-2000:-200),iterations = 10^5,shape=adap.al,output="${HOME}/rates.jld2",);
```

If you have a system with few resources, it is possible [to download pre-computed TGP and DGN data rates](https://imkt.uab.cat/files/inputs/rates.jld2). Please go to the next section to continue the inference using the pre-computed rates.