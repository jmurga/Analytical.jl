using Pkg

const ENV["R_HOME"]="*"

Pkg.add("Conda")

using Conda

Conda.add("r-base",channel="conda-forge")
Conda.add(["r-locfit","r-ggplot2","r-data.table","r-r.utils"],channel="conda-forge")

Pkg.add("RCall")