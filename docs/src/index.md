# Analytical

Analytical approximation to $\alpha_{x}$ accounting for linkage. We explore the impact of linkage and background selection at positive selected alleles sites. The package  solve anylitical approximations for different genetic scenarios in order to estimate the strenght and rate of adaptation. 

When empircal values of polymorphim and divergence are given, they will be used to discern their expected correspoding values modeled under any Distribution of Fitness Effect (*DFE*) and background selection values (*B*). 

Our goal is to subset summary statistics given a set of empirical values for any genetic scenario that would be used as prior distributions in *ABC* algorithms.


```@docs
Analytical.fixNeut
Analytical.fixNegB
Analytical.pFix
Analytical.fixPosSim
Analytical.DiscSFSNeutDown
Analytical.poissonFixation
Analytical.poissonPolymorphism2
Analytical.alphaByFrequencies
Analytical.phiReduction
adap
```

# Installation

To install our module we highly recommend to use [LTS official Julia binaries](https://julialang.org/downloads/). If is your first time using Julia, you can easily export the Julia bin through ```export PATH="/path/to/directory/julia-1.0.5/bin:$PATH"``` in your shell. Since we use *scipy* to solve equations, you must install it on your Python. Once Julia and *scipy* are installed you just need to run at your shell:

```bash
julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
```

Or from Pkg REPL (by pressing `]` at Julia interpreter):

```julia
add PackageSpec(path="https://github.com/jmurga/Analytical.jl")
```

In addition we provide a Docker image based on Debian including Julia and Jupyter notebook. You can access to Debian system or just to Jupyter pulling the image from dockerhub. Remember to link the folder `/analysis` with any folder at your home to save the results:

```bash
# Pull the image
docker pull jmurga/mktest
# Run docker bash interactive session linking to some local volume to export data. 
docker run -i -t -v ${HOME}/folderPath:/analysis/folder  jmurga/mktest
# Run only jupyter notebook from docker image. Change the port if 8888 is already used
docker run -i -t -v ${HOME}/folderPath:/analysis/folder -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```



# Dependencies

# References
