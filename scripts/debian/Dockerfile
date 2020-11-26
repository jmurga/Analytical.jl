FROM debian:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel r-base && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

ENV LD_LIBRARY_PATH /opt/miniconda/lib:$LD_LIBRARY_PATH

RUN conda update -n base -c defaults conda
RUN conda init bash
RUN conda install -c conda-forge jupyterlab numpy scipy mpmath pandas -y

RUN R -e 'install.packages(c("abc","ggplot2","data.table","dplyr"), repos = "https://cloud.r-project.org",dependecies=TRUE)'

RUN wget --quiet https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.1-linux-x86_64.tar.gz -P /root
RUN  tar -zxvf /root/julia-1.5.1-linux-x86_64.tar.gz -C /root

ADD data/mk_with_positions_BGS.txt /analysis/mk_with_positions_BGS.txt
# ADD scripts/precompile_mktest.jl /root/precompile_mktest.jl

RUN /root/julia-1.5.1/bin/julia -e 'using Pkg;Pkg.add(["IJulia","Parameters", "NLsolve", "SpecialFunctions", "Distributions", "Roots", "ArbNumerics", "StatsBase", "LsqFit", "PoissonRandom", "SparseArrays","OrderedCollections","GZip","DataFrames","Parsers","CSV","PackageCompiler","RCall"])'

RUN /root/julia-1.5.1/bin/julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
# RUN /root/julia-1.5.1/bin/julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:Analytical, sysimage_path="/root/mktest.so", precompile_execution_file="/root/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /root/ABCreg
RUN cd /root/ABCreg/src && make

ENV PATH /root/julia-1.5.1/bin:$PATH
ENV PATH /root/ABCreg/src:$PATH

WORKDIR /analysis
