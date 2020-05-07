# Home

Analytical estimation of $\alpha_{x}$ accounting for linkage. We explore the impact of linkage and background selection at positive selected alleles sites. The package  solve anylitical approximations for different genetic scenarios in order to estimate the strenght and rate of adaptation. 

When empircal values of polymorphim and divergence are given, they will be used to discern their expected correspoding values modeled under any Distribution of Fitness Effect ($DFE$) and background selection values ($B$). 

Our goal is to subset summary statistics given a set of empirical values for any genetic scenario that would be used as prior distributions in $ABC$ algorithms.


```@docs
Analytical.fixNeut
Analytical.fixNegB
Analytical.pFix
Analytical.fixPosSim
Analytical.DiscSFSNeutDown
```

