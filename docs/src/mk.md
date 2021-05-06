# MK approaches
We included other heuristic MK approaches in our module. All the functions use the formated SFS and divergence data described at [data section](data.md).

## Standard MKT
The standard McDonald & Kreitman test (MKT) (McDonald and Kreitman, 1991) was developed to be applied to protein coding sequences, combining both divergence (D) and polymorphism (P) sites, and categorizing mutations as synonymous (P0, D0) and non-synonymous (Pi, Di). 

If all mutations are either strongly deleterious or neutral, then $Di/D0$ is expected to roughly equal $Pi/P0$. In contrast, if positive selection is operating in the region, adaptive mutations rapidly reach fixation and contribute more to divergence than polymorphism compared to neutral mutations, and then $Di/D0 > Pi/P0$. Assuming that adaptive mutations contribute little to polymorphism but substantially to divergence, the proportion of non-synonymous substitutions is inferred following Smith and Eyre-Walker (2002)

$\alpha = 1 - (\frac{P_i}{P_0}\cdot\frac{D_0}{D_i}$ 


## Fay, Waycoff and Wu MKT
Fay et al. (2001) proposed an approach that removes all polymorphisms segregating at a frequency below a given threshold (usually 5%–15%). Although there is no consensus about what this threshold should be used, J. Charlesworth & Eyre-Walker (2008) demonstrated that  estimates are robust using a frequency threshold of 15%, below which most slightly deleterious polymorphisms are found and removed. The estimates are reasonably accurate only when the rate of adaptive evolution is high and the Distribution of Fitness Effects (DFE) of deleterious mutations is leptokurtic (J. Charlesworth & Eyre-Walker, 2008).

## imputed MKT (in preparation)

## Asymptotic MKT
Proposed by Messer and Petrov (2013). This extension is robust to the presence of selective sweeps (genetic hitchhiking) and the segregation of slightly deleterious polymorphisms substitutions (BGS). In this approach, the authors defined $\alpha$ as a function that depends on the SFS of alleles. Therefore, $\alpha$ is estimated in different frequency intervals ($x$), and these values are then adjusted to an exponential function. An exponential fit is suitable as the non-synonymous allele frequency is expected to decay exponentially over the respective levels of synonymous polymorphisms (Messer & Petrov, 2013).

$\alpha_{fit(x)} = a+b \cdot e^{-cx}$