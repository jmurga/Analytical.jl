# ABC-MK

ABC-MK is an analytical approximation to $\alpha_{(x)}$. We have explored the impact of linkage and background selection at positive selected alleles sites. The package solves analytical approximations for different genetic scenarios to estimate the strength and rate of adaptation.

Our approach directly estimates $\alpha_{(x)}$ and several statistics ($B$, $\alpha_W$, $\alpha_S$) associated with random DFE. In conjunction, the associated values to these DFE can be used as summary statistics at ABC methods. Therefore, our method can estimate the rate and strength of adaption in models and non-models organisms.

## Docker installation
We highly recommend using the Docker image to execute the software. The Docker image is based on Debian and includes all the software needed to run the pipeline. You can access to Debian system or Jupyter pulling the image from [Docker Hub](https://hub.docker.com/r/jmurga/abcmk). Remember to link the folder /analysis with any folder at your ```${HOME}``` directory to save the results:


```bash
# Pull the image
docker pull jmurga/abcmk
# Run docker linking some local volume to export data
docker run -i -t -v analysis/:/analysis/folder jmurga/abcmk
# Run jupyter notebook from docker image. Change the port if 8888 is already used
docker run -i -t -v analysis/:/analysis/folder -p 8888:8888 jmurga/abcmk /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```

To use our command-line interface, just run

```bash
docker run -i -t -v ${HOME}/test/:/analysis/test/ jmurga/abcmk /analysis/abcmk_cli.jl
```

## Singularity installation
We have created a Singularity container to use the software in HPC systems. We have tested the software at HPC working with Slurm and HTCondor scheduler

```singularity
singularity pull --arch amd64 library://jmurga/default/abcmk:latest
```

We found a [bug](https://github.com/JuliaParallel/ClusterManagers.jl/issues/164) regarding Singularity, ClusterManagers.jl and Slurm in our HPC tests. Please, consider to install the packages manually if your HPC works with Slurm. We provided a [Julia script](https://github.com/jmurga/Analytical.jl/blob/master/scripts/julia_dependencies.jl) to easily install all the required packages. Just run it before to execute our [Command-Line Interface](cli.jl). We provide [Slurm](https://github.com/jmurga/Analytical.jl/blob/master/scripts/abcmkSlurm.sh) and [HTCondor](https://github.com/jmurga/Analytical.jl/blob/master/scripts/abcmkHtcondor.sub) examples showing the  estimation.

## Scratch installation
To install our module from scratch, we highly recommend using [LTS official Julia binaries](https://julialang.org/downloads/). Once you have installed Julia in your system, consider to install some important dependencies to automatize your ABC-MK pipelines. We have prepared a file to install them. Please download [this file](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/julia_dependencies.jl) and execute the following command

```bash
curl -o julia_dependencies.jl https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/julia_dependencies.jl
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
# Linux debian-based installation
sudo apt install libgsl-dev libz-dev build-essential git
# MacOS installation
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install GSL zlib git
```

```bash
git clone https://github.com/molpopgen/ABCreg.git ABCreg
cd ABCreg/src && make
```

### R
We have used R to estimate the Maximum-A-Posteriori (MAP) from posterior distributions following ABCreg examples. We linked Julia and R internally. The module contains functions to perform the estimations without quit the Julia session.

```bash
# Linux debian-based installation
sudo apt install r-base
# MacOS installation
brew install r
```

If you are going to perform MAP estimates and plot using our module, be sure you have installed R and the following packages: ggplot2 and data.table, locfit. 

```R
R -e "install.packages(c('ggplot2','data.table','locfit'))"
```

## References
- Uricchio, L.H., Petrov, D.A. & Enard, D. Exploiting selection at linked sites to infer the rate and strength of adaptation. Nat Ecol Evol 3, 977–984 (2019). [https://doi.org/10.1038/s41559-019-0890-6](https://doi.org/10.1038/s41559-019-0890-6)
- Philipp W. Messer, Dmitri A. Petrov. Frequent adaptation and the McDonald–Kreitman test. Proceedings of the National Academy of Sciences May 2013, 110 (21) 8615-8620. [https://doi.org/10.1073/pnas.1220835110](https://doi.org/10.1073/pnas.1220835110)
- Nordborg, M., Charlesworth, B., & Charlesworth, D. (1996). The effect of recombination on background selection. Genetical Research, 67(2), 159-174. [https://doi.org/10.1017/S0016672300033619](https://doi.org/10.1017/S0016672300033619)
- R R Hudson and N L Kaplan. Deleterious background selection with recombination. [Genetics December 1, 1995 vol. 141 no. 4 1605-1617.](https://www.genetics.org/content/141/4/1605)
- Linkage and the limits to natural selection. N H Barton. [Genetics June 1, 1995 vol. 140 no. 2 821-841](https://www.genetics.org/content/140/2/821)
- Thornton, K.R. Automating approximate Bayesian computation by local linear regression. BMC Genet 10, 35 (2009). [https://doi.org/10.1186/1471-2156-10-35](https://doi.org/10.1186/1471-2156-10-35)
