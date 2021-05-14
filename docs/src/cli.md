# Command-Line Interface

We develop a Command-Line Interface (CLI) in case you want to avoid Julia interpreter. You can easly download [abcmk_cli.jl](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/abcmk_cli.jl). The CLI have different functions to perform the whole pipeline as explained in [Infering the rate and strength of adaptation](empirical.md) section

```bash
julia abcmk_cli.jl  
```

```bash
See --help of each command for usages  
  rates  
  parseTgpData  
  parseDgnData  
  summaries  
  abcInference  
  plotMap  
```

To reproduce the examples you can follow the steps described at [Empirical estimation](https://jmurga.github.io/Analytical.jl/dev/empirical/#Computational-pipeline-1) section

## Estimating rates
To perform the rate estimations you can use the function ```rates``` at the CLI. The function works simirarly to the function [```Analytical.rates```](@ref). You can check the argument at [Rates](https://jmurga.github.io/Analytical.jl/dev/rates/#Estimating-fixation-and-polymorphic-rates-considering-generalized-model-of-selection-and-linkage-1) sectino or using the macro ```@doc Analytical.rates``` in the Julia interpreter.

```julia
julia abcmk_cli.jl rates --help
```

```
julia abcmk_cli.jl rates --samples 661 --gamNeg -1000,-200 --gL 1,10 --gH 400,1000 --rho 0.001 --theta 0.001 --solutions 1000000 --output rates.jld2 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000  --nthreads 7

Function to solve fixation and polymorphic rates analitically. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.

If rho and/or theta are set to nothing, the function will input random values given the range 0.0005:0.0005:0.01. Otherwise you can fix the values. If gL is set to nothing, the function will not account the role of the weakly selected alleles in the estimation.

The function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.

Please check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/Analytical.jl/dev/cli/ to check model

Optional Arguments:
  --ne: Int64 (default: 1000)
  --samples: Int64 (default: 500)
  --gamNeg: String (default: -1000,-200)
  --gL: String (default: 5,10)
  --gH: String (default: 400,1000)
  --dac: String (default: 2,4,5,10,20,50,200,500,700)
  --shape: Float64 (default: 0.184)
  --rho: String (default: nothing)
  --theta: String (default: nothing)
  --solutions: Int64 (default: 100000)
  --output: String (default: /home/jmurga/rates.jld2)
  --scheduler: String (default: local)
  --nthreads: Int64 (default: 1)

```

If you are going to perform the estimation in a HPC, please set the variable ```scheduler``` using the name of the HPC task manager. By default the value is set to ```local```

```bash
time julia abcmk_cli.jl rates --samples 661 --gamNeg -2000,-200 --gL 1,10 --gH 200,2000 --rho 0.001 --theta 0.001 --solutions 100000 --output /home/jmurga/rates.jld2 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000 --nthreads 7 --scheduler local
```

## Parse data into new folder
To estimate summary statistics, you need to provide empirical SFS and divergence files. As explained in section [data](data.md), you can directly parse TGP or DGN data using our module. Nonetheless, you can input any other SFS and divergence file.

The function ```parseData``` will create a folder containing TGP or DGN rawdata, parsed SFS and divergence files. Otherwise you can specify an exisiting folder.

```bash
julia abcmk_cli.jl parseData --help
```

```bash
julia abcmk_cli.jl parseData --analysisFolder analysis/ --geneList analysis/dnaVipsList.txt

Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. 

The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics.	

Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli/

Optional Arguments:
  --analysisFolder: String (default: <folder>)
  --dataset: String (default: tgp)
  --geneList: String (default: false)
  --bins: String (default: false)
```

```bash
julia abcmk_cli.jl parseData --analysisFolder analysis/
```

Remember you can use the argument ```geneList``` to subset genes from TGP or DGN data using a list of Ensembl or Flybase ID. Please check [Multiple dataset](https://jmurga.github.io/Analytical.jl/dev/multiple/) to get more info.

## Estimate summary statistics

To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md).


```bash
julia abcmk_cli summaries --help
```

```bash
julia abcmk_cli.jl summaries --analysisFolder analysis/ --rates analysis/rates.jld2 --samples 661 --dac 2,4,5,10,20,50,200,661,925 --summstatSize 1000000 --nthreads 7


Estimate summary statistics from analytical rates. You must provide a path containing the parsed SFS and divergence file.

The function returns files containing bootstrapped datasets (alphas.txt) and summary statistics (summstat.txt)

Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli

Optional Arguments:
  --analysisFolder: String (default: <folder>)
  --rates: String (default: rates.jld2)
  --ne: Int64 (default: 1000)
  --samples: Int64 (default: 500)
  --dac: String (default: 2,4,5,10,20,50,200,500,700)
  --summstatSize: Int64 (default: 100000)
  --replicas: Int64 (default: 100)
  --bootstrap: String (default: true)
  --scheduler: String (default: local)
  --nthreads: Int64 (default: 1)
```

The function will output observed data bootstraped (*alphas.txt*) and summary statistics (*summaries.txt*) in the analysisFolder. These file will be used at ABC inference to generate posterior distributions.

```bash
julia abcmk_cli.jl summaries --analysisFolder /home/jmurga/tgpData/ --rates /home/jmurga/rates.jld2 --samples 661 --replicas 100 --summstatSize 100000 --dac 2,4,5,10,20,50,200,661,925 --nthreads 7
```

## Perform ABC inference
At this point, you have a folder containing summary statistics and observed data to perform ABC inference. As explained in our [home page](index.md), we performed the ABC inference using [ABCreg](https://github.com/molpopgen/ABCreg). However, you can used other ABC software to perform the inference.

We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. If you are going to use ABCreg to make inferences from our software directly, please [cite the publication](https://doi.org/10.1186/1471-2156-10-35). Remember you need to install ABCreg before continue. Please check [home page](index.md) to install ABCreg.

It is possible to perform the inference through Julia. The number of parameters to infer will always be 5: $\alpha_w$, $\alpha_s$, $\alpha$, $\gamma$ and $\beta$. We set the tolerance value 
 The function will output one file per bootstrapped replicas containing the posteriors distributions. The posterior distributions contains five columns corresponding to:
 - $\alpha_w$: Contribution of weak selecction to $\alpha$
 - $\alpha_s$: Contribution of strong selecction to $\alpha$
 - $\alpha$: Adaptation rate
 - $\gamma$: Negative selection coefficient


```bash
julia abcmk_cli.jl abcInference --analysisFolder /home/jmurga/tgpData/ --replicas 100 --P 5 --S 9 --tol 0.001 --ABCreg /home/jmurga/ABCreg/src/reg --parallel true --nthreads 7
```

## Estimate Maximum-A-Posteriori and plot using R. Using julia expression, cannot input into *abcmk_cli.jl* (in development)

We used R to estimate the Maximum-A-Posteriori (MAP) from posterior distributions following ABCreg examples. We linked Julia and R internally. The module contains functions to perform the estimations without quit the Julia session.

If you will perform MAP estimates and plot using our module, be sure you have installed R and the following packages: ggplot2 and data.table, locfit. 


```bash
julia -e 'using Analytical, RCall, GZip, DataFrames, CSV;Analytical.sourcePlotMapR(script="/home/jmurga/tgpData/script.jl"); Analytical.plotMap(analysisFolder="/home/jmurga/tgpData");'
```