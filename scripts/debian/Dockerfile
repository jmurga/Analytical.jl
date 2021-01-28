FROM julia

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel r-base && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("abc","ggplot2","data.table","dplyr"), repos = "https://cloud.r-project.org",dependecies=TRUE)'

ADD data/mk_with_positions_BGS.txt /analysis/mk_with_positions_BGS.txt

RUN julia -e 'using Pkg;Pkg.add(["IJulia","PackageCompiler" ,"Parameters", "NLsolve", "SpecialFunctions", "Distributions", "Roots", "ArbNumerics", "StatsBase", "LsqFit", "PoissonRandom", "SparseArrays", "Distributed", "CSV", "SharedArrays", "JLD2","DataFrames","GZip","Parsers","OrderedCollections","FastaIO","RCall"])'

ADD scripts/precompile_mktest.jl /root/precompile_mktest.jl

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
RUN julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:Analytical, sysimage_path="/root/mktest.so", precompile_execution_file="/root/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /root/ABCreg
RUN cd /root/ABCreg/src && make

# ENV PATH /root/julia-1.5.1/bin:$PATH
ENV PATH /root/ABCreg/src:$PATH

WORKDIR /analysis

# CMD ["julia"]
CMD ["julia","--sysimage","/root/mktest.so"]