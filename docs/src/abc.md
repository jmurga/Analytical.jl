# *ABC* inference from empirical data

## Prior distributions from analytical estimations

In this example we provide a solution to replicate results at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). We will simulate $10^6$ summary statistics from random *DFE* to use as prior distribution in *ABCreg*. In this case we will need a set of empirical observed values in order to subset the summary statistics.

We need to set the model accounting for the sampling value. The *SFS* is expected to be in raw frequencies. If the model is not properly set up, the *SFS* will not be correctly parsed. In our case, we are going to set up a model with default parameters only to parse the *SFS* and convolute the observed frequencies with a binomial distribution.

```julia
using Analytical
adap = Analytical.parameters(N=1000,n=661)
```

Once the model account for the number of samples we can open the files. The function `Analytical.parseSfs` will return polymorphic and divergent counts and SFS accounting for the whole spectrum: `collect(1:adap.nn)/adap.nn`. In addition an output file will be created contained the observed values to input in *ABCreg*.

```julia
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

pol, sfs, d = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/testData",sfsColumns=[3,5],divColumns=[6,7],bins=100)
```

The module include a function to solve *N* times different genetic scenarios. We solve the analytical approximation taking into account random and independent values to draw *DFE* and $\alpha_{(x)}$. Each parameter combination are replicated to 5% frequency bins background selection values (saved at `adap.bRange`).
```julia
# Execute one to compile the function
Analytical.summaryStats(param=adap,alpha=0.4,divergence=d,sfs=sfs,bins=100,iterations=1);
# Make your estimations
df = Analytical.summaryStats(param=adap,alpha=0.4,divergence=d,sfs=sfs,bins=100,iterations=10^5);
```

To parallelize the process we created a thread pool inside [`summaryStats`](@ref) using the *Distributed* package. To parallel the process you only need to define the available process and add our model to each thread.

```julia
using Distributed
nthreads=4; addprocs(nthreads)
# Load the module in all the threads
@everywhere using Analytical, DataFrames, CSV
adap = Analytical.parameters(N=1000,n=661)
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

pol, sfs, d = Analytical.parseSfs(param=adap,data=files,output="/home/jmurga/testData",sfsColumns=[3,5],divColumns=[6,7],bins=100)

# Execute one to compile the function
Analytical.summaryStats(param=adap,alpha=0.4,divergence=d,sfs=sfs,bins=100,iterations=1);
# Make your estimations
df = Analytical.summaryStats(param=adap,alpha=0.4,divergence=d,sfs=sfs,bins=100,iterations=10^5);
CSV.write("/home/jmurga/prior", DataFrame(df), delim='\t',header=false);

```

Alternatively, if you are going to use the command line script, please make the threads available when executing julia
```bash
julia --procs 8 script.jl --arg1 --arg2 --arg3
```

## *ABC* inference
Generic ABC methods proceed by three main steps: (1) first sampling parameter values from prior distributions, (2) next simulating a model and calculating informative summary statistics, and lastly (3) comparing the simulated summary statistics to observed data. The parameter values that produce summary statistics that best match the observed data form an approximate posterior distribution. We link Julia with *ABCreg*. It will output one file per line in data. The files contain the posterior distributions. We return the posterior distributions, mean and quantiles.

```julia
posteriors, results  = Analytical.ABCreg(data="/home/jmurga/data.tsv",prior="/home/jmurga/prior.tsv", nparams=27, nsummaries=24, outputPath="/home/jmurga/", outputPrefix="outPaper", tolerance=0.001, regressionMode="T",regPath="/home/jmurga/ABCreg/src/reg")
```

You can easily plot the posterior distributions using Julia or just input the files at your favorite plot software.

```julia
using Plots, Plots.PlotMeasures, StatsPlots, LaTeXStrings

function plotPosterior(data,file,imgSize)

	Plots.gr()
	Plots.theme(:wong2)

	p1 = StatsPlots.density(data],
							legend    = :topright,
							fill      = (0, 0.3),
							xlabel    = L"\alpha",
							label     = [L"\alpha_S" L"\alpha_W" L"\alpha"],
							ylabel    = "Posterior density",
							lw        = 0.5,
							fmt       = :svg,
							size      = imgSize
						)
	Plots.savefig(file)

	return p1
end

p = plotPosterior(posterior[1],"/home/jmurga/posterior1.svg",(800,600))
```

[![NBViewer](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jmurga/Analytical.jl/blob/master/scripts/analyticalAlphaAndPriors.ipynb)
