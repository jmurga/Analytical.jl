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

Declare a variable containing some basic information about your model. We used a sample size of 661 to perform later analysis over TGP data. The selected Derived Alleles Counts (DAC) will be used to compute summary statistics and perform ABC inference. It is possible to subset any of the selected DAC values when computing summary statistics.

```julia
adap = Analytical.parameters(n=661,dac=[1,2,4,5,10,20,50,100,200,400,500,661,925,1000])
```

Note that ```adap``` contains about mutation rate, recombination rates, DFE, BGS and probabilities of fixations. To check all the arguments you can access to the function documentation using ```@doc Analytical.parameter```

```julia
@doc Analytical.parameter

  Mutable structure containing the variables required to solve the analytical approach. All the functions are solve using the internal values of the structure. For this reason, adap is the
  only exported variable. Adap should be change before the perform the analytical approach, in other case, $\alpha_{(x)}$ will be solve with the default values.

  Parameters
  ≡≡≡≡≡≡≡≡≡≡≡≡

    •  gamNeg::Int64: Selection coefficient for deleterious alleles

    •  gL::Int64: Selection coefficient for weakly benefical alleles

    •  gH::Int64: Selection coefficient for strongly benefical alleles

    •  alLow::Float64: Proportion of α due to weak selection

    •  alTot::Float64: α

    •  thetaF::Float64: Mutation rate defining BGS strength

    •  thetaMidNeutral::Float64: Mutation rate on coding region

    •  al::Float64: DFE shape parameter

    •  be::Float64: DFE scale parameter

    •  B::Float64: BGS strength

    •  bRange::Array{Float64,1}: BGS values to simulate

    •  pposL::Float64: Fixation probabily of weakly beneficial alleles

    •  pposH::Float64: Fixation probabily of strongly beneficial alleles

    •  N::Int64: Population size

    •  n::Int64: Sample size

    •  Lf::Int64: Flanking region length

    •  rho::Float64: Recombination rate

    •  TE::Float64
```

We will estimate the empirical adaptation rate using TGP data using the estimated DFE parameters at Boyko et al (2008). Nonetheless, shape and scale DFE parameter are flexible in our model.

```julia
adap.al = 0.184
adap.be = abs(0.184/-457)
```
Before to automatize the fixation and polimorphic rates estimation, you must to convolute the binomial distribution to obtain the downsampled SFS

```julia
convolutedSamples = Analytical.binomialDict()
Analytical.binomOp!(adap,convolutedSamples.bn)
```

Note the ```Analytical.binomOp!``` make inplace estimation at ```convolutedSamples```given the BGS range defined at ```adap.bRange```

```julia
@docs Analytical.binomOp!

  binomOp(param)

  Binomial convolution to sample the allele frequencies probabilites depeding on background selection values, and sample size.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  param::parameters

    •  convolutedSamples::binomialDict

  Returns
  ≡≡≡≡≡≡≡≡≡

    •  Array{Float64,2}: convoluted SFS given for each B value defined in the model. Results saved at param.bn.

```


Now the variable ```adap``` contains sample size, DAC and DFE information. The function ```Analytical.rates``` will perform the analytical estimation of *N* independent models regarding DFE, BGS, mutation rate, and recombination rate. In the following example, we used the function ```Analytical.rates``` to input the prior distributions. The function will randomize the input values to solve *N* independent estimation regarding our model. 

```julia
@doc Analytical.rates
```

```julia
@time df = Analytical.rates(param = adap,convolutedSamples=convolutedSamples,gH=collect(200:2000),gL=collect(1:10),gamNeg=collect(-2000:-200),iterations = 10^5,shape=adap.al,output="analysis/rates.jld2");
```

The function will create a HDF5 file containing the solved models, fixation rates, polymorphic rates, and the selected DAC. This information will be used later to estimate summary statistics.


Note that [```Analytical.rates```](@ref) is the most resource and time-consuming function. In our case, the function will estimate 10^5 independent models. Each model solves the estimation for all possible BGS values. We used BGS values from 0.1 to 0.999 in 5% increments. In total, the example will produce 3.7 million estimates. We have used a hierarchical data structure (HDF5) to facilitate model parameters and rates storage.

The following example took about 1.5 hours to execute on the hardware described at section [Infering the rate and strength of adaptation](empirical.md)

If you have a system with few resources, it is possible [to download pre-computed TGP and DGN data rates](https://imkt.uab.cat/files/inputs/rates.jld2). Please go to the next section to continue the inference using the pre-computed rates.