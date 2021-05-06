# *ABC* inference from empirical data

At this point, you have a folder containing summary statistics and observed data to perform ABC inference. As explained in our [home page](index.md), we performed the ABC inference using [ABCreg](https://github.com/molpopgen/ABCreg). However, you can used other ABC software to perform the inference.

We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform *ABC* inference. If you are going to use *ABCreg* to make inferences from our software directly, please [cite the publication](https://doi.org/10.1186/1471-2156-10-35). Remember you need to install ABCreg before continue. Please check [home page](index.md) to install ABCreg.

It is possible to perform the inference through Julia. The number of parameters to infer will always be 5: $\alpha_w$, $\alpha_s$, $\alpha$, $\gamma$ and $\beta$

```julia
Analytical.ABCreg(analysisFolder="/home/jmurga/tgpData/",replicas=100,P=5,S=size(adap.dac,1),tol=0.002,workers=1,abcreg="/home/jmurga/ABCreg/src/reg",parallel=false);
```

If you have installed GNU-parallel in your system it is possible to parallelize the inference

```julia
Analytical.ABCreg(analysisFolder="/home/jmurga/tgpData/",replicas=100,P=5,S=size(adap.dac,1),tol=0.002,workers=7,abcreg="/home/jmurga/ABCreg/src/reg",parallel=true);
```

We used R to estimate the Maximum-A-Posteriori (MAP) from posterior distributions following ABCreg examples. We linked Julia and R internally. The module contains functions to perform the estimations without quit the Julia session.

If you will perform MAP estimates and plot using our module, be sure you have installed R and the following packages: ggplot2 and data.table, locfit. 

```julia
Analytical.sourcePlotMapR(script="/home/jmurga/tgpData/script.jl")
tgpmap = Analytical.plotMap(analysisFolder="/home/jmurga/tgpData/");
describe(tgpmap)
```

```
 Row │ variable  mean          min           median        max          nmissing  eltype   
     │ Symbol    Float64       Float64       Float64       Float64      Int64     DataType 
─────┼─────────────────────────────────────────────────────────────────────────────────────
   1 │ aw           0.108927     0.010857       0.108796      0.19754          0  Float64
   2 │ as           0.0506607   -0.00750128     0.0514826     0.134143         0  Float64
   3 │ a            0.152842     0.0962341      0.149083      0.233131         0  Float64
   4 │ gamNeg    1184.81       512.458       1277.47       1903.12             0  Float64
   5 │ shape        0.142934     0.128369       0.14189       0.167394         0  Float64
```

![image](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/docs/src/figure2.svg =500x500)