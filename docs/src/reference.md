## Model parameters

*adap* is the only variable exported from *Analytical* module. It is a Mutable structure contaning the variables required to solve the analytical approach. Any value can be easly changed. Remember *adap* should be change before the execution, in other case, $\alpha_{(x)}$ will be solve with the default values. To change all the values at once, you can use [`Analytical.parameters`](@ref) in order to set specific models in a new variable.

```@docs
Analytical.parameters
Analytical.binomialDict
Analytical.Br
Analytical.setThetaF!
Analytical.setPpos!
Analytical.binomOp!
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
Analytical.DiscSFSSelPosDown
Analytical.DiscSFSSelNegDown
Analytical.cumulativeSfs
```

## Summary statistics
```@docs
Analytical.analyticalAlpha
Analytical.poissonFixation
Analytical.poissonPolymorphism
Analytical.sampledAlpha
Analytical.alphaByFrequencies
Analytical.bgsIter
Analytical.summaryStats
Analytical.gettingRates
Analytical.summaryStatsFromRates
Analytical.ratesToSummaries
```

## Inference tools
```@docs
Analytical.parseSfs
Analytical.ABCreg
```

## Raw data tools
```@docs
Analytical.uSfsFromFasta
```
