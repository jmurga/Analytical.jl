# ABC-MK

ABC-MK is an analytical approximation to $\alpha_{(x)}$. We have explored the impact of linkage and background selection at positive selected alleles sites. The package solves analytical approximations for different genetic scenarios to estimate the strength and rate of adaptation.

Our approach directly estimates $\alpha_{(x)}$ and several statistics ($B$, $\alpha_W$, $\alpha_S$) associated with random DFE. In conjunction, the associated values to these DFE can be used as summary statistics at ABC methods. Therefore, our method can estimate the rate and strength of adaption in models and non-models organisms.

## Docker installation
We highly recommend using the Docker image to execute the software. The Docker image is based on Debian and includes all the software needed to run the pipeline. You can access to Debian system or Jupyter pulling the image from [Docker Hub](https://hub.docker.com/r/jmurga/abcmk). Remember to link the folder /analysis with any folder at your ```${HOME}``` directory to save the results:


```bash
# Pull the image
docker pull jmurga/abcmk
# Run docker linking some local volume to export data
docker run -i -t -v ${HOME}/<folderPath>:/analysis/folder jmurga/abcmk
# Run jupyter notebook from docker image. Change the port if 8888 is already used
docker run -i -t -v ${HOME}/<folderPath>:/analysis/folder -p 8888:8888 jmurga/abcmk /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```

To use our command-line interface, just run

```bash
docker run -i -t -v ${HOME}/test/:/analysis/test/ jmurga/abcmk /analysis/abcmk_cli.jl
```

## Singularity installation
We have created a Singularity container to use the software in HPC systems

```singularity
singularity pull --arch amd64 library://jmurga/default/abcmk:latest
```

## Scratch installation
To install our module from scratch, we highly recommend using [LTS official Julia binaries](https://julialang.org/downloads/)

```bash
JULIA_VERSION=1.6.1
curl -o ${HOME}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-${JULIA_VERSION}-linux-x86_64.tar.gz
tar -zxf ${HOME}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C ${HOME}
``` 

If this is your first time using Julia, please remember to export Julia binaries to your path. In this way, you will execute Julia using the command ```julia```

```bash
echo 'export PATH="${HOME}/julia-${JULIA_VERSION}/bin:$PATH"' >> ${HOME}/.bashrc
source ${HOME}/.bashrc
```

Consider to install some important dependencies to automatize your ABC-MK pipelines. We have prepared a file to install them. Please download [this file](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/julia_dependencies.jl) and execute the following command

```bash
wget https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/julia_dependencies.jl
julia julia_dependencies.jl
```

You can easly install our Julia package executing

```bash
julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
```

Or from Pkg REPL (by pressing `]` at Julia interpreter):

```julia
add https://github.com/jmurga/Analytical.jl
```

### ABCreg
We have linked [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. Nonetheless others ABC softwares could be used ([abc (R package)](https://doi.org/10.1111/j.2041-210X.2011.00179.x), [ABCToolBox](https://doi.org/10.1186/1471-2105-11-116), etc). If you are going to use ABCreg to directly make inference from our software please [cite the publication](https://doi.org/10.1186/1471-2156-10-35) and compile it in your system. Anyway, once you get the summary statistic files you can use any other ABC software.

ABCreg needs *GSL* and *libz* to work. Please install both libraries before compile the software:

```bash 
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

### R
We have used R to estimate the Maximum-A-Posteriori (MAP) from posterior distributions following ABCreg examples. We linked Julia and R internally. The module contains functions to perform the estimations without quit the Julia session.

If you are going to perform MAP estimates and plot using our module, be sure you have installed R and the following packages: ggplot2 and data.table, locfit. 

```R
R -e "install.packages(c('ggplot2','data.table','locfit'))"
```

## Dependencies
All Julia dependencies are installed within the package. You don't need to install them manually.

### Mandatory dependencies to solve the analytical equations
- [`Roots`](https://github.com/JuliaMath/Roots.jl) - root finding.
- [`NLsolve`](https://github.com/JuliaStats/Distributions.jl) - non-linear systems of equations solver.
- [`Distributions`](https://github.com/JuliaStats/Distributions.jl) - probability distributions.
- [`SpecialFunctions`](https://github.com/JuliaMath/SpecialFunctions.jl) - special mathematical functions in Julia.
- [`Quadmath`](https://github.com/JuliaMath/Quadmath.jl) - multiprecision numerical computing.
- [`PoissonRandom`](https://github.com/SciML/PoissonRandom.jl) - Poisson random number generator.
- [`Parameters`](https://github.com/mauro3/Parameters.jl) - custom keyword constructor.

### The following dependencies are required to use all the funcionalities (parse SFS, plots, parse multi-Fasta, etc.)
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
