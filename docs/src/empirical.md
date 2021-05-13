# Infering the rate and strength of adaptation

We develop an extension of the analytical approximation presented in Uricchio, Petrov, and Enard (2019). In the previous paper, analytical calculations were used to explore the effect of BGS on weakly beneficial alleles, but the estimation procedure employed was based on computationally intensive simulations. While inference procedures based on forward simulations allow for the inclusion of complex demographic histories, they often require high-performance computing resources and can be prohibitively slow for some models. Here we extend Uricchio, Petrov and Enard (2019) analytical approximations to develop a simple and computationally efficient ABC-based inference procedure. Our method accounts for the DFE of deleterious and beneficial alleles and incomplete recombination between selected genomic elements. 


To perform the empirical estimation of $\alpha_{(x)}$ we followed a generic ABC algorithm (Beaumont et al 2002). ABC proceeds by first sampling parameter values from prior distributions, the next simulating model using these parameter values, calculating informative summary statistics, and comparing the simulated summary statistics to the observed data. Considering the standard ABC scheme, we (1) considered $N$ random combinations from prior distributions; (2) solved $N$ independent models to generate informative summary statistics; (3) subset the parameter values producing summary statistics that best match the observed data from an approximate posterior distribution. Additionally, a linear model can be imposed to correct the non-0 distance between the simulated and observed summary statistics.

## Extending analytical estimations
We extended the analytical calculations solving $N$ independents models. We automatize the analytical estimations to input prior distributions. 

The models will be solved using random combinations from prior distributions. Fixation, polymorphic rates, and model information will be used to estimate informative summary statistics needed to perform ABC inference. Please see section [rates](rates.md) to check how to input prior distributions.

## Summary statistics
\cite{uricchio_exploiting_2019} used the analytical theory to demonstrate the effect of weakly beneficial alleles at $\alpha_{(x)}$. To do that, they solved $\alpha_{(x)}$ through the fixation and polymorphism rates since the locus length ($L$) and the time branch ($T$) estimations at first order $\alpha_{(x)}$ estimation can be canceled

$\mathbb{E}[\alpha_{(x)}] \approx 1 - \frac{LT(\mathbb{E}[d_S])}{LT(\mathbb{E}[d_N])} \frac{LT(\mathbb{E}[p_{S(x)}])}{LT(\mathbb{E}[p_{N(x)}])} \approx 1 - \frac{\mathbb{E}[d_S]}{\mathbb{E}[d_N]} \frac{\mathbb{E}[p_{S(x)}]}{\mathbb{E}[p_{N(x)}]}$

where $\frac{LT}{LT}$ with a constant mutation rate tend to $1$. 

Here, we extended their analytical calculations to generate the summary statistics required at ABC approaches. Thereby we avoided expensive forward simulations. 

Considering fixed values of $T$, $L$, $\mu$, and fixation rates it is possible to solve analytical $\alpha_{(x)}$ performing the multiplication. However, this requires not only branch length estimations but also explicitly locus length selection which highly increases the order of estimations to solve. To avoid branch length estimations and locus length selection, we follow the previously described assumptions: (1) empirically observed fixations should be proportional to $T$, $L$ and $\mu$; (2) the mutational process follows a Poisson distribution . Based on that premises, we used a Poisson-sampling process to estimate the expected number of fixations. We corrected the analytical expected rates by the empirical observations as the rate of success $\lambda$ on the Poisson distribution

$\mathbb{E}[D] = X \in Poisson\left(\lambda = D_{observed} \times (\mathbb{E}[d_+]+\mathbb{E}[d_-]+\mathbb{E}[d_0])\right)$

To estimate $\alpha_{(x)}$ we draw the expected values considering each fixation category. We sampled the values in two categories according to the relative rates. Following the same procedure:

$\mathbb{E}[D_S] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_0]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

$\mathbb{E}[D_N] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_+] + \mathbb{E}[d_-]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

Our model takes into account that both sampling variance and process variance should affect the number of variable alleles that we sample at any particular allele frequency. The process variance arises from the random mutation-fixation process along the branch. To incorporate that variance we made one sample per frequency-bin given the SFS. Considering the random process at each frequency we made credible intervals for future ABC estimations. Otherwise, we would falsely find higher confidence in our parameter estimates. We draw the expected polymorphic values similarly to fixations, considering the SFS and the expected rates.

$\mathbb{E}[P] = \sum_{x=0}^{x=1} X \in Poisson\left(\lambda = SFS_{(x)\ observed} \times (\mathbb{E}[p_{+(x)}]+\mathbb{E}[p_{-(x)}]+\mathbb{E}[p_{0(x)}])\right)$

We modified the expected polymorphic rates and the observed SFS considering each frequency $x$ as the cumulative sum above $x$. This quantities have the same asymptote but are less affected by changing sample size. Therefore, we scaled better the analysis to sample size since most common alleles at frequencies $x$ have very few polymorphic sites in large samples. This way, we finally transformed $\alpha_{(x)}$ to be depending on the previous frequency bin value, which reflects each frequency category over the expected asymptotic shape. Because of that, $\alpha_{(x)}$ analysis is more robust even in low polymorphic populations.

Consequently, for each analytical combination, it is possible to sample the expected number of fixations and polymorphism to perform the $\alpha_{(x)}$ estimations. Like in \cite{uricchio_exploiting_2019} we input $\alpha_{(x)}$ as summary statistic at generic ABC algorithm. We exploited summary statistics selection that are informative for estimating $\alpha_{(x)}$ values.

## Parameters inference
We used the expected proportion of weakly and strong fixations to estimate $\alpha_W$, $\alpha_S$ and $\alpha$. We followed the same procedure as the above section to subset the expected number of fixation taking into account weakly, strong and deleterious fixation rates categories. Therefore, we modified \hyperref[eqn 19]{eqn 19} according to their relative fixation rates.

$\mathbb{E}[D_W] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_W]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

$\mathbb{E}[D_S] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_S]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

In addition, for each model, we retrieve negative selection coefficients and shape parameter to perform ABC inference.

# Computational pipeline
The following sections describe a pipeline to estimate the empirical adaptation rate ($\alpha$) using whole-genome protein, Viral Interacting Proteins and Non-Viral Interacting Protein data from the TGP. This analysis is similar to the one performed at Uricchio et. al (2019). We divided the pipeline into the following step:

 - Analytical estimations
 - Empirical data parse
 - Summary statistic estimation
 - ABC inference

Please note that we use a relative path called *analysis/* to execute the whole pipeline. This path will be relative to the where Julia or the Command-Line Interface is executed.

The software is prepared to parallelize each pipeline step using Julia [Distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) computing. Distributing the process into threads has a cost in RAM. Please make some tests in your machine before executing expensive models. Nonetheless, distributing the estimation into threads decreases the pipeline execution time. It is almost mandatory to parallelize at least the rate estimations.

The following examples, as well as the analysis described at REF, were tested using a laptop with the following hardware:
- Intel i7-7700HQ (8) @ 3.800GHz 
- 16GB RAM DDR4 2400MHz
