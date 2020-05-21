# Analytical estimation.
### Solving \alpha_(x)
Our method is based in the analytical solution of $\alpha_{(x)}$ given a genetic scenario. The method could be extended to *DFE* and background selection values in order to get summary statistics that can be used as prior distrubtions at *ABC* methods. In this example we show how asympotic $\alpha$ is affected by linkage and background selection.

*adap* is the only variable exported from *Analytical* module. It is a Mutable structure contaning the variables required to solve the analytical approach. Any value can be easly changed. Remember *adap* should be change before the execution, in other case, $\alpha_{(x)}$ will be solve with the default values. To change all the values at once, you can use [`Analytical.changeParameters`](@ref) in order to set specific models. Please take into account some package used here could be changed. 

Load the modules
```julia
using Analytical, Plots
```

Set the model parameters
```julia
Analytical.changeParameters(gam_neg=-83,gL=10,gH=500,alLow=0.2,alTot=0.2,theta_f=1e-3,theta_mid_neutral=1e-3,al=0.184,be=0.000402,B=0.999,bRange=append!(collect(0.2:0.05:0.95),0.999),pposL=0.001,pposH=0.0,N=500,n=25,Lf=10^6,rho=0.001,TE=5.0,convoluteBinomial=true)
```

Solving the model. Here we solve $\alpha$ generally using the expected rates. We are not considering any especific mutation process over a locus and branch time.

```julia
Analytical.set_theta_f()
theta_f = adap.theta_f
adap.B = 0.999
Analytical.set_theta_f()
Analytical.setPpos()
adap.theta_f = theta_f
adap.B= B

x,y = Analytical.analyticalAlpha(gammaL=adap.gL,gammaH=adap.H,pposL=adap.pposL,pposH=adap.pposH)
```

Plotting the results. $x$ contains $\alpha_{(x)}$ accounting for weakly beneficial alleles. $y$ contains the true value of $\alpha_{(x)}$, not accounting for weakly beneficial alleles.

```julia
Plots.PlotMeasures
Plots.gr()
Plots.theme(:bright)

Plots.plot(collect(1:size(x,1)),hcat(x,y),
    legend = :bottomright,
    seriestype = :line,
    xlabel = "Derived Alleles Counts",
    ylabel = "α",
    label = ["All alleles" "Neutral + deleterious"],
    markershapes= [:circle :circle],
    lw = 1,
    xscale = :log,
    fmt = :svg,
    bottom_margin=10 * Plots.PlotMeasures.mm,
    left_margin=10 * Plots.PlotMeasures.mm,
    size=(700,500)
)
```

![image](https://raw.githubusercontent.com/jmurga/Analytical.jl/aa980e3de977ba5e9c3f536fb4d5459ef8035221/docs/src/fig1.svg)


### Empirical case
In this example we provide a solution to replicate results at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). We will simulate $10^6$ summary statistics from random *DFE* to use as prior distribution in *ABCreg*. In this case we will need a set of empirical observed values in order to subset the summary statistics.

We need to set the model accounting for the sampling value. The SFS is expected to be in raw frequencies. If the model is not properly set up, the SFS will not be correctly parse. In our case, we are going to set up a model with default parameters only to parse the SFS and convolute the observed frequencies with a binomial distribution.

```julia
Analytical.changeParameters(N=1000,n=661,convoluteBinomial=true)
```

Once the model account for the number of samples we can open the files. The function `Analytical.parseSfs` will return polymorphic and divergent counts and SFS accounting for the whole spectrum: `collect(1:adap.nn)/adap.nn`. In addition an output file will be created contained the observed values to input in *ABCreg*.

```julia
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

empiricalValues = Analytical.parseSfs(data=files,output="testData.tsv",sfsColumns=[3,5],divColumns=[6,7])
```