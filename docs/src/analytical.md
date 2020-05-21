# Analytical estimation.
### Solving \alpha_(x)
Our method is based in the analytical solution of $\alpha_{(x)}$ given a genetic scenario. The method could be extended to *DFE* and background selection values in order to get summary statistics that can be used as prior distrubtions at *ABC* methods. In this example we show how asympotic $\alpha}$ is affected by linkage and background selection.

*adap* is the only variable exported from *Analytical* module. It is a Mutable structure contaning the variables required to solve the analytical approach. Any value can be easly changed. Remember *adap* should be change before the execution, in other case, $\alpha_{(x)}$ will be solve with the default values. To change all the values at once, you can use [`Analytical.changeParameters`](@ref) in order to set specific models.

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

Ploting the resutlts. $x$ contains $\alpha_{(x)}$ accounting for weakly beneficial alleles. $y$ contains the true value of $\alpha_{(x)}$, not accounting for weakly beneficial alleles.

```julia
Plots.PlotsMeasures
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


![image](/docs/fig1.svg)
### Empirical case
In this example we provide a solution to replicate results at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). Please taking into account some package used here could be changed. 

```julia
using Analytical, BenchmarkTools, Plots
```