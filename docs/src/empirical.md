# Infering the rate and strength of adaptation

We develop an extension of the analytical approximations presented in Uricchio, Petrov, and Enard (2019). In the previous paper, analytical calculations were used explore the effect BGS on weakly beneficial alleles, but the estimation procedure employed was based on computationally intensive simulations. While inference procedures based on forward simulations allow for the inclusion of complex demographic histories, they often require high performance computing resources and can be prohibitively slow for some models. Here we extend analytical approximations of Uricchio, Petrov and Enard (2019) and develop a simple and computationally efficient ABC-based inference procedure that accounts for the DFE of deleterious and beneficial alleles and incomplete recombination between selected genomic elements. 

The next sections describe the whole pipeline to automatize analytical estimations, empircal data parse, summary statistic estimation and ABC inference.

The software is prepared to parallelize each pipeline step using Julia [Distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) computing. Distributing the process into threads has a cost in RAM memory, please make some tests in your machine before to execute expensives models. Nonetheless, distributing the estimation into threads highly decrease the pipeline execution time. It is almost mandatory to parallelize at least the rates estimations.

The following examples as well as the analysis described at REF were tested using a laptop with the following hardware:
- Intel i7-7700HQ (8) @ 3.800GHz 
- 16GB RAM DDR4 2400MHz
