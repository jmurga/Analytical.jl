# Analytical estimation
## Solving $\alpha_{(x)}$
Our method is based in the analytical solution of $\alpha_{(x)}$ given a genetic scenario. The method could be extended to any *DFE* and background selection values in order to get summary statistics that can be used at *ABC* methods. In this example we show how the asympotic $\alpha$ is affected by linkage and background selection.

Before start you need to set up a variable of type *Analytical.parameters*. It is a Mutable structure containing the variables required to solve the analytical approach. Any value can be easily changed. Remember that the variable should be defined before the execution, in other case, $\alpha_{(x)}$ will be solve with default values. To change all the values at once, just check the variables at the struct [`Analytical.parameters`](@ref) in order to set specific models.

#### Load the modules
```julia
using Analytical
```

#### Set the model parameters and convolute the binomial distribution.
We set a global adaptation rate of 0.4 with a contribution of 0.2 regarding weak adaptation. The process is modeled using selection coefficient of 500,10 for strongly and weakly beneficial alleles. We modeled the DFE for deleterious alleles using the values provided at [Boyko et al. (2008)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000083). 

```julia
adap = Analytical.parameters(n=661,alTot=0.4,alLow=0.2,gH=500,gamNeg=-457,al=0.184,be = 0.184/457,B=0.999)
```

Before solve the model, we compute the convoluted binomial distribution to subset the expected SFS

```julia
convolutedSamples = Analytical.binomialDict()
Analytical.binomOp!(adap,convolutedSamples.bn)
```

#### Solve the model
Here we solve $\alpha_{(x)}$ generally, using the expected rates. We are not considering any specific mutation process over a locus and branch time.

```julia
x,y = Analytical.analyticalAlpha(param = adap, convolutedSamples = convolutedSamples)
```

Internally the function first set the mutation rate regarding the BGS strenght and the probability of fixations given the genetic scenario. Then it estimate the SFS and fixations rates. Please check [`analyticalAlpha`](@ref) to understand the process.

#### Plotting the results.
$x$ contains $\alpha_{(x)}$ accounting for weakly beneficial alleles. $y$ contains the value of $\alpha_{(x)}$, not accounting for weakly beneficial alleles. In this example, we do not model BGS (check *B* parameter in *Analytical.parameters* function above)

In Julia you can easily use R ggplot packages using the module [RCall](). Please check you have installed R in your system. Nonetheless you can plot using Julia Plots module.

Install RCall just in case you want to plot using R

```julia
using Pkg
Pkg.add("RCall")
``` 

```julia
using RCall

R"""
	library(ggplot2)
	library(reshape2)
	df = data.frame(f=seq(1,length($x)),x=$x,y=$y)
	d  = melt(df,id.vars='f')
	p = ggplot(df,aes(x=f,y=value,color=variable)) + geom_line() + geom_point() + scale_colour_manual(values=c('#30504f', '#e2bd9a'),labels = c("Nuetral + deleterious alleles", "All alleles")) + theme_bw()
	p
"""
```

![image](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/docs/src/fig1.svg)
