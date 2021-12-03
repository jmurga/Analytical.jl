# Command-Line Interface

We develop a Command-Line Interface (CLI) in case you want to avoid Julia interpreter. You can easly download [abcmk_cli.jl](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/abcmk_cli.jl). The CLI have different functions to perform the whole pipeline as explained in [Infering the rate and strength of adaptation](empirical.md) section

```bash
julia abcmk_cli.jl  
```

```bash
See --help of each command for usages  
  rates  
  parse_data
  summaries  
  inference
```

To reproduce the examples you can follow the steps described at [Empirical estimation](https://jmurga.github.io/Analytical.jl/dev/empirical/#Computational-pipeline-1) section

## Estimating rates
To perform the rate estimations you can use the function ```rates``` at the CLI. The function works simirarly to the function [```Analytical.rates```](@ref). You can check the argument at [Rates](https://jmurga.github.io/Analytical.jl/dev/rates/#Estimating-fixation-and-polymorphic-rates-considering-generalized-model-of-selection-and-linkage-1) sectino or using the macro ```@doc Analytical.rates``` in the Julia interpreter.

```julia
julia abcmk_cli.jl rates --help

usage: abcmk_cli.jl rates [--pop_size POP_SIZE]
                        --sample_size SAMPLE_SIZE --dac DAC
                        --gam_neg GAM_NEG
                        --positive_strong POSITIVE_STRONG
                        --positive_weak POSITIVE_WEAK [--shape SHAPE]
                        [--rho RHO] [--theta THETA]
                        [--solutions SOLUTIONS] [--output OUTPUT]
                        [--scheduler SCHEDULER] [--nthreads NTHREADS]
                        [-h]

optional arguments:
  --pop_size POP_SIZE   Population size (type: Int64, default: 1000)
  --sample_size SAMPLE_SIZE
                        Sample size (type: Int64)
  --dac DAC             Derived Allele Count
  --gam_neg GAM_NEG     Selection coefficient for deleterious alleles
  --positive_strong POSITIVE_STRONG
                        Selection coefficient for strongly beneficial
                        alleles
  --positive_weak POSITIVE_WEAK
                        Selection coefficient for weakly beneficial
                        alleles
  --shape SHAPE         Shape value modeling Gamma distribution for
                        deleterious alleles (type: Float64, default:
                        0.184)
  --rho RHO             Recombination rate (type: Float64, default:
                        0.001)
  --theta THETA         Mutation rate on coding locus (type: Float64,
                        default: 0.001)
  --solutions SOLUTIONS
                        Mutation rate on coding locus (type: Int64,
                        default: 100000)
  --output OUTPUT       Output file (default: "rates.jld2")
  --scheduler SCHEDULER
                        Select scheduler manager (default: "local")
  --nthreads NTHREADS   Select number of threads to parallelize (type:
                        Int64, default: 1)
  -h, --help            show this help message and exit

```

If you are going to perform the estimation in a HPC, please set the variable ```scheduler``` using the name of the HPC task manager. By default the value is set to ```local```

```bash
time julia abcmk_cli.jl rates --sample_size 661 --dac 1,2,4,5,10,20,50,100,200,500,661,925 --gam_neg -2000:-200 --positive_strong 200:1000 --positive_weak 1:100 --shape 0.184 --output rates.jld2 --solutions 100000 --scheduler local --nthreads 7
```

## Parse data into new folder
To estimate summary statistics, you need to provide empirical SFS and divergence files. As explained in section [data](data.md), you can directly parse TGP or DGN data using our module. Nonetheless, you can input any other SFS and divergence file.

The function ```parse_data``` will create a folder containing TGP or DGN rawdata, parsed SFS and divergence files. Otherwise you can specify an exisiting folder.

```bash
julia abcmk_cli.jl parse_data --help
Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. 

The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics. 

Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli/
```

```bash
julia abcmk_cli.jl parse_data --gene_list analysis/dnaVipsList.txt analysis/

Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. 

The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics.	

Please check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli/

positional arguments:
  folder                a positional argument

optional arguments:
  -d, --dataset DATASET
                        a positional argument (default: "tgp")
  -g, --gene_list GENE_LIST a positional argument (type: Union{Bool, String}, default: false)
  -b, --bins BINS       a positional argument (type: Union{Bool,Int64}, default: false)
  -h, --help            show this help message and exit


```

```bash
julia abcmk_cli.jl parse_data analysis/
```

Remember you can use the argument ```geneList``` to subset genes from TGP or DGN data using a list of Ensembl or Flybase ID. Please check [Multiple dataset](https://jmurga.github.io/Analytical.jl/dev/multiple/) to get more info.

## Estimate summary statistics

To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md).


```bash
julia abcmk_cli summaries --help

Estimate summary statistics from analytical rates. You must provide a path containing the parsed SFS and divergence file.

The function returns files containing bootstrapped datasets (alphas.txt) and summary statistics (summstat.txt)

Check the documentation to get more info https://jmurga.github.io/Analytical.jl/dev/cli

usage: abcmk_cli.jl summaries --rates RATES --sample_size SAMPLE_SIZE
                        --dac DAC --summstat_size SUMMSTAT_SIZE
                        [--bootstrap BOOTSTRAP] [--replicas REPLICAS]
                        [--nthreads NTHREADS] [--scheduler SCHEDULER]
                        [-h] folder

positional arguments:
  folder                Folder path containing SFS and divergence files to run the analysis
optional arguments:
  --rates RATES         H5 file containing precomputed rates
  --sample_size SAMPLE_SIZE Sample size (type: Int64)
  --dac DAC             Derived allele count
  --summstat_size SUMMSTAT_SIZE Define number of summary estatistic to perform ABC (type: Int64)
  --bootstrap BOOTSTRAP Allow bootstrap following polyDFE manual (type: Bool, default: false)
  --replicas REPLICAS   Number of bootstrap replicas (type: Int64, default: 1)
  --nthreads NTHREADS   Select scheduler manager (type: Int64, default: 1)
  --scheduler SCHEDULER Select scheduler manager (default: "local")
  -h, --help            show this help message and exit

```

The function will output observed data bootstraped (*alphas.txt*) and summary statistics (*summaries.txt*) in the analysisFolder. These file will be used at ABC inference to generate posterior distributions.

```bash
julia abcmk_cli.jl summaries --rates rates.jld2 --sample_size 661 --dac 2,4,5,10,20,50,100,661,925 --summstat_size 100000 --bootstrap true --replicas 100 --nthreads 10 --scheduler andromeda analysis/
```

## Perform ABC inference
At this point, you have a folder containing summary statistics and observed data to perform ABC inference. As explained in our [home page](index.md), we performed the ABC inference using [ABCreg](https://github.com/molpopgen/ABCreg). However, you can used other ABC software to perform the inference.

We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. If you are going to use ABCreg to make inferences from our software directly, please [cite the publication](https://doi.org/10.1186/1471-2156-10-35). Remember you need to install ABCreg before continue. Please check [home page](index.md) to install ABCreg.

It is possible to perform the inference through Julia. The function will output one file per bootstrapped replicas containing the posteriors distributions. We set the tolerance value to record 1000 values for the regression.  The posterior distributions contains five columns corresponding to:

 - α weak: Contribution of weak selecction to $\alpha$
 - α strong: Contribution of strong selecction to $\alpha$
 - α weak: Adaptation rate
 - γ: Negative selection coefficient
 - β: Negative selection coefficient

```bash
julia abcmk_cli.jl inference --help

positional arguments:
  folder                Folder path containing SFS and divergence files to run the analysis

optional arguments:
  --S S                 Define number of summary estatistic to perform
                        ABC (type: Int64)
  --tol TOL             Tolerance (type: Float64)
  --abcreg ABCREG       ABCreg path
  --nthreads NTHREADS   Select scheduler manager (type: Int64, default: 1)
  --scheduler SCHEDULER Select scheduler manager (default: "local")
  -h, --help            show this help message and exit



```

```bash
julia abcmk_cli.jl inference --S 9 --tol 0.01 --abcreg /home/jmurga/ABCreg/src/reg --nthreads 7 analysis/
```

## Estimate Maximum-A-Posteriori and plot using R. 

Using julia expression, cannot input into *abcmk_cli.jl* (in development)

We used R to estimate the Maximum-A-Posteriori (MAP) from posterior distributions following ABCreg examples. We linked Julia and R internally. The module contains functions to perform the estimations without quit the Julia session.

If you will perform MAP estimates and plot using our module, be sure you have installed R and the following packages: ggplot2 and data.table, locfit. 


```bash
julia -e 'using Analytical, RCall, GZip, DataFrames, CSV;Analytical.source_plot_map_R(script="analysis/script.jl"); Analytical.plot_map(analysisFolder="analysis/");'
``` 