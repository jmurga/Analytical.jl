# *adap* variable
*adap* is the only variable exported from *Analytical* module. It is a Mutable structure contaning the variables required to solve the 4 analytical approach. Any value could be easily change reassing. *adap* should be change before the execution, in other case, $\alpha_{(x)$ will be solve with the default values. To change all the values at once, you can use [Analytical.changeParameters](@ref Analytical.changeParameters) in order to set a model.

```@docs
Analytical.adap
```
# Analytical estimation
## Fixations
```@docs
Analytical.fixNeut
Analytical.fixNegB
Analytical.pFix
Analytical.fixPosSim
```
## Polymorphism
```@docs
Analytical.DiscSFSNeutDown
```

## Summary statistics
```@docs
Analytical.poissonFixation
Analytical.poissonPolymorphism
Analytical.alphaByFrequencies
```