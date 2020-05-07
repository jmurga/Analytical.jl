# Analytical

[![Build Status](https://travis-ci.com/jmurga/Analytical.jl.svg?branch=master)](https://travis-ci.com/jmurga/Analytical.jl) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jmurga.github.io/Analytical.jl/stable)  

Extended ABC-MK calculations accounting for background selection and weak adaptation based on Julia language.

Included in this repository are the scripts used to simulate and infer parameters from ABC software.

To install the module we highly recommend to use [LTS official Julia binaries](https://julialang.org/downloads/). You can easily export the Julia bin through ```export PATH="/path/to/directory/julia-1.0.5/bin:$PATH"```. In addition, since the package use *scipy* functions, we recommend executing Julia activating the following [conda enviroment](https://github.com/jmurga/Analytical.jl/tree/master/scripts/abc-mk.yml) (or you can just just install *scipy* in your default python). Once Julia is installed just run:
```bash
julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
```

We provide a Docker image based on Debian including Julia and all packages needed to run the estimations. It also includes jupyterlab and several data science packages. You can run the whole Debian image or just the jupyterlab instance pulling the image from dockerhub:

```bash
# Pull the image
docker pull jmurga/mktest
# Run docker bash interactive session linking to some local volume to export data
docker run -i -t -p 8888:8888 -v ${HOME}/<anyData>:/data/<anyData> -e HOSTID=$(id -u) jmurga/mktest
# Run only jupyter notebook from docker image. Change the port if 8888 is already used
docker run -i -t -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```
