# Analytical estimation
## Solving $\alpha_{(x)}$
Our method is based in the analytical solution of $\alpha_{(x)}$ given a genetic scenario. The method could be extended to any *DFE* and background selection values in order to get summary statistics that can be used as prior distrubtions at *ABC* methods. In this example we show how asympotic $\alpha$ is affected by linkage and background selection.

Before start you need to set up a variable of type *Analytical.parameters*. It is a Mutable structure containing the variables required to solve the analytical approach. Any value can be easily changed. Remember that the variable should be defined before the execution, in other case, $\alpha_{(x)}$ will be solve with default values. To change all the values at once, just check the variables at the struct [`Analytical.parameters`](@ref) in order to set specific models. Please take into account some package used here could be changed.

#### Load the modules
```julia
using Analytical, Plots
```

#### Set the model parameters and convolute the binomial distribution.
```julia
adap = Analytical.parameters(gam_neg=-83,
                            gL=10,
                            gH=500,
                            alLow=0.2,
                            alTot=0.2,
                            theta_f=1e-3,
                            theta_mid_neutral=1e-3,
                            al=0.184,
                            be=0.000402,
                            B=0.999,
                            bRange=append!(collect(0.2:0.05:0.95),
                            0.999),
                            pposL=0.001,
                            pposH=0.0,
                            N=500,
                            n=25,
                            Lf=10^6,
                            rho=0.001,
                            TE=5.0,
                        )
Analytical.binomOp!(adap);
```

#### Solve the model
Here we solve $\alpha$ generally using the expected rates. We are not considering any specific mutation process over a locus and branch time.

```julia
x,y = Analytical.analyticalAlpha(gammaL=adap.gL,gammaH=adap.H,pposL=adap.pposL,pposH=adap.pposH)
```

Internally the function first set the mutation rate and the probability of fixations given the genetic scenario. Then it estimate the SFS and fixations rates. Please check [`AnalyticalAlpha`](@ref) to understand the process.

#### Plotting the results.
$x$ contains $\alpha_{(x)}$ accounting for weakly beneficial alleles. $y$ contains the true value of $\alpha_{(x)}$, not accounting for weakly beneficial alleles.

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

![image](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/docs/src/fig1.svg)
