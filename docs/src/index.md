# Analytical

Analytical approximation to $\alpha_{x}$. We explore the impact of linkage and background selection at positive selected alleles sites. The package solves anylitical approximations for different genetic scenarios in order to estimate the strenght and rate of adaptation. 

When empircal values of polymorphim and divergence are given, they will be used to discern their expected correspoding values modeled under any Distribution of Fitness Effect (*DFE*) and background selection values (*B*). 

Our goal is to subset summary statistics given empirical observed values that would be used as prior distributions in *ABC* algorithms. Please check

# Installation

To install our module we highly recommend to use [LTS official Julia binaries](https://julialang.org/downloads/). If is your first time using Julia, you can easily export the Julia bin through ```export PATH="/path/to/directory/julia-1.v.v/bin:$PATH"``` in your shell. Since we use *scipy* to solve equations, the package depends on PyCall.

```bash
julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
```

Or from Pkg REPL (by pressing `]` at Julia interpreter):

```julia
add https://github.com/jmurga/Analytical.jl
```

## Scipy installation
You can install *scipy* on your default Python or install it through Julia Conda:

```julia
julia -e 'using Pkg;Pkg.add("PyCall");using PyCall;pyimport_conda("scipy.optimize", "scipy")'
```

If you cannot install properly *scipy* through Julia Conda try the following:

- Set an empty Python enviroment and re-build PyCall: `ENV["PYTHON"]="";  Pkg.build("PyCall")`
- Re-launch Julia and install the scipy.optimize module: `using PyCall;pyimport_conda("scipy.optimize", "scipy)`

## Docker
We provide a Docker image based on Debian including Julia and Jupyter notebook. You can access to Debian system or just to Jupyter pulling the image from dockerhub. Remember to link the folder `/analysis` with any folder at your home to save the results:

```bash
# Pull the image
docker pull jmurga/mktest
# Run docker bash interactive session linking to some local volume to export data. 
docker run -i -t -v ${HOME}/folderPath:/analysis/folder  jmurga/mktest
# Run only jupyter notebook from docker image. Change the port if 8888 is already used
docker run -i -t -v ${HOME}/folderPath:/analysis/folder -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```

# Dependencies
All the dependecies are installed within the package. You don't need to install manually. If you experiment any problem contact us or try to reinstall *Pycall* and *scipy*.

### Mandatory dependencies to solve the analytical equations
- [`Roots`](https://github.com/JuliaMath/Roots.jl) - root finding.
- [`Distributions`](https://github.com/JuliaStats/Distributions.jl) - probability distributions.
- [`PyCall`](https://github.com/JuliaPy/PyCall.jl) - directly call and fully interoperate with Python.
- [`SpecialFunctions`](https://github.com/JuliaMath/SpecialFunctions.jl) - special mathematical functions in Julia.
- [`Parameters`](https://github.com/mauro3/Parameters.jl) - custom keyword constructor.


### The following dependencies are required to use all the funcionalities (parse SFS, plots, etc.)
- [`CSV`](https://github.com/JuliaNLSolvers/Optim.jl)
- [`Parsers`](https://github.com/JuliaStats/Distributions.jl)
- [`StatsBase`](https://github.com/JuliaStats/Distances.jl)
- [`DataFrames`](https://github.com/JuliaStats/Distances.jl)
- [`GZip`](https://github.com/JuliaIO/GZip.jl)

### ABC
We link [ABCreg](https://github.com/molpopgen/ABCreg) with julia in order to perform ABC estimations. If you are going to use ABCreg to do inference please [cite the publication](https://doi.org/10.1186/1471-2156-10-35) and compile it in your system. Anyway, once you get the priors distributions you can use any other ABC software.

```bash
git clone https://github.com/molpopgen/ABCreg.git
cd ABCreg/src && make
```

# References
