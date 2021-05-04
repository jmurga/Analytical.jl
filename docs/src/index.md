# ABC-MK

ABC-MK is an analytical approximation to $\alpha_{(x)}$. We explore the impact of linkage and background selection at positive selected alleles sites. The package solves analytical approximations for different genetic scenarios in order to estimate the strength and rate of adaptation.

Our approach estimates directly $\alpha_{(x)}$ and several statistics ($B$, $\alpha_W$, $\alpha_S$) associated to random *DFE*.  In conjunction, the associated values to these *DFE* can be used as priors distributions at *ABC* methods. If we subset enough parameters, we will consider any frequency spectrum and fixations under generalized models of selection, demography, and linkage associated with the empirical population and sample size. Therefore, our method can estimate rate and strength of adaption in models and non-models organisms, for which previous *DFE* and demography are unknown.

## Installation
### Docker
We highly recommend to use the Docker image to execute the software. The Docker image is based on Debian and include all the software needed to run the pipeline. You can access to Debian system or Jupyter pulling the image from dockerhub. Remember to link the folder `/analysis` with any folder at your ```${HOME}``` directory to save the results:

```bash
# Pull the image
docker pull jmurga/abcmk
# Run docker bash interactive session linking to some local volume to export data.
docker run -i -t -v ${HOME}/<folderPath>:/analysis/folder  jmurga/mktest
# Run only jupyter notebook from docker image. Change the port if 8888 is already used
docker run -i -t -v ${HOME}/<folderPath>:/analysis/folder -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```

You can directly perform the analysis without using the interactive command-line running the script [abcmk_cli.jl](https://github.com/jmurga/Analytical.jl/blob/master/scripts/abcmk_cli.jl):
```bash
julia abcmk_cli.jl --help
```

or using Docker (remember to link some folder to output your results!)
```bash
docker run -i -t -v ${HOME}/folderPath:/analysis/folder jmurga/abcmk abcmk_cli.jl --help
```

### Scratch installation
To install our module from scratch we highly recommend to use [LTS official Julia binaries](https://julialang.org/downloads/)

```
JULIA_VERSION=1.6.1
curl -o ${HOME}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-${JULIA_VERSION}-linux-x86_64.tar.gz
tar -zxf ${HOME}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C ${HOME}
``` 

If is your first time using Julia please remember to export Julia binaries to your path. In this way you will execute Julia using the command ```julia```

```bash
echo 'export PATH="${HOME}/julia-${JULIA_VERSION}/bin:$PATH"' >> ${HOME}/.bashrc
source ${HOME}/.bashrc
```

You can easly install the Julia package executing

```bash
julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
```

Or from Pkg REPL (by pressing `]` at Julia interpreter):

```julia
add https://github.com/jmurga/Analytical.jl
```

## ABC
We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia in order to perform *ABC* inference, although another *ABC* software could be used ([abc (R package)](https://doi.org/10.1111/j.2041-210X.2011.00179.x), [ABCToolBox](https://doi.org/10.1186/1471-2105-11-116), etc). If you are going to use *ABCreg* to directly make inference from our software please [cite the publication](https://doi.org/10.1186/1471-2156-10-35) and compile it in your system. Anyway, once you get the summary statistic files you can use any other *ABC* software.

ABCreg needs *GSL* and *libz* to work. Please install both libraries before compile the software:

```bash 
# Installing make and git just in case they are not available in your system
sudo apt install libgsl-dev libz-dev build-essential git
```

```bash
git clone https://github.com/molpopgen/ABCreg.git ${HOME}/ABCreg
cd ${HOME}/ABCreg/src && make
```

Once ABCreg is installed you can add the binary to your PATH:

```bash 
echo 'export PATH="${HOME}/ABCreg/src/:$PATH"' >> ${HOME}/.bashrc
source ${HOME}/.bashrc
```

## Dependencies
All the dependencies are installed within the package. You don't need to install manually.

#### Mandatory dependencies to solve the analytical equations
- [`Roots`](https://github.com/JuliaMath/Roots.jl) - root finding.
- [`NLsolve`](https://github.com/JuliaStats/Distributions.jl) - non-linear systems of equations solver.
- [`Distributions`](https://github.com/JuliaStats/Distributions.jl) - probability distributions.
- [`SpecialFunctions`](https://github.com/JuliaMath/SpecialFunctions.jl) - special mathematical functions in Julia.
- [`Quadmath`](https://github.com/JuliaMath/Quadmath.jl) - multiprecision numerical computing.
- [`PoissonRandom`](https://github.com/SciML/PoissonRandom.jl) - Poisson random number generator.
- [`Parameters`](https://github.com/mauro3/Parameters.jl) - custom keyword constructor.


#### The following dependencies are required to use all the funcionalities (parse SFS, plots, parse multi-Fasta, etc.)
- [`CSV`](https://github.com/JuliaNLSolvers/Optim.jl)
- [`Parsers`](https://github.com/JuliaStats/Distributions.jl)
- [`StatsBase`](https://github.com/JuliaStats/Distances.jl)
- [`DataFrames`](https://github.com/JuliaStats/Distances.jl)
- [`GZip`](https://github.com/JuliaIO/GZip.jl)
- [`OrderedCollections`](https://github.com/JuliaCollections/OrderedCollections.jl)
- [`Plots`](https://github.com/JuliaPlots/Plots.jl)
- [`FastaIO`](https://github.com/carlobaldassi/FastaIO.jl)

## References
- Uricchio, L.H., Petrov, D.A. & Enard, D. Exploiting selection at linked sites to infer the rate and strength of adaptation. Nat Ecol Evol 3, 977–984 (2019). [https://doi.org/10.1038/s41559-019-0890-6](https://doi.org/10.1038/s41559-019-0890-6)
- Philipp W. Messer, Dmitri A. Petrov. Frequent adaptation and the McDonald–Kreitman test. Proceedings of the National Academy of Sciences May 2013, 110 (21) 8615-8620. [https://doi.org/10.1073/pnas.1220835110](https://doi.org/10.1073/pnas.1220835110)
- Nordborg, M., Charlesworth, B., & Charlesworth, D. (1996). The effect of recombination on background selection. Genetical Research, 67(2), 159-174. [https://doi.org/10.1017/S0016672300033619](https://doi.org/10.1017/S0016672300033619)
- R R Hudson and N L Kaplan. Deleterious background selection with recombination. [Genetics December 1, 1995 vol. 141 no. 4 1605-1617.](https://www.genetics.org/content/141/4/1605)
- Linkage and the limits to natural selection. N H Barton. [Genetics June 1, 1995 vol. 140 no. 2 821-841](https://www.genetics.org/content/140/2/821)
- Thornton, K.R. Automating approximate Bayesian computation by local linear regression. BMC Genet 10, 35 (2009). [https://doi.org/10.1186/1471-2156-10-35](https://doi.org/10.1186/1471-2156-10-35)
