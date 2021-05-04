# *ABC* inference from empirical data

## Prior distributions from analytical estimations

In this example we provide a solution to replicate results at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). We will simulate $3.7 \cdot 10^6$ summary statistics from random *DFE* to use as prior distribution in *ABCreg*. In this case we will need a set of empirical observed values in order to subset the summary statistics.

We need to set the model accounting for the sampling value. The *SFS* is expected to be in raw frequencies. If the model is not properly set up, the *SFS* will not be correctly parsed. In our case, we are going to set up a model with default parameters only to parse the *SFS* and convolute the observed frequencies with a binomial distribution.

```julia
using Analytical
adap = Analytical.parameters(n=661)
```

Once the model account for the number of samples we can open the files. The function `Analytical.parseSfs` will return polymorphic and divergent counts and SFS accounting for the whole spectrum: `collect(1:adap.nn)/adap.nn`. In addition an output file will be created contained the observed values to input in *ABCreg*.

```julia

```

The module include a function to solve *N* times different genetic scenarios. We solve the analytical approximation taking into account random and independent values to draw *DFE* and $\alpha_{(x)}$. Each parameter combination are replicated to 5% frequency bins background selection values (saved at `adap.bRange`).

```julia
# Execute one to compile the function
```

To parallelize the process we created a thread pool inside [`summaryStats`](@ref) using the *Distributed* package. To parallel the process you only need to define the available process and add our model to each thread.

```julia
# Load Distributed package and add threads

# Load the pacakgein all the threads
```

Alternatively, if you are going to use the command line script, please make the threads available when executing julia
```bash
```

## *ABC* inference
Generic ABC methods proceed by three main steps: (1) first sampling parameter values from prior distributions, (2) next simulating a model and calculating informative summary statistics, and lastly (3) comparing the simulated summary statistics to observed data. The parameter values that produce summary statistics that best match the observed data form an approximate posterior distribution. We link Julia with *ABCreg*. It will output one file per line in data. The files contain the posterior distributions. We return the posterior distributions, mean and quantiles.


You can easily plot the posterior distributions using Julia or just input the files at your favorite plot software.

[![NBViewer](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.jupyter.org/github/jmurga/Analytical.jl/blob/master/scripts/analyticalAlphaAndPriors.ipynb)
