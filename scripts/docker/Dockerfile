FROM julia

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel r-base && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

ENV PATH /opt/conda/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

RUN conda install -c conda-forge jupyter -y

RUN conda init bash 

RUN R -e 'install.packages(c("abc","ggplot2","data.table","dplyr","IRkernel"), repos = "https://cloud.r-project.org",dependecies=TRUE)'
RUN R -e 'IRkernel::installspec()'

ADD data/mk_with_positions_BGS.txt /analysis/mk_with_positions_BGS.txt

RUN julia -e 'using Pkg;Pkg.add(["Fire","IJulia","PackageCompiler" ,"Parameters", "NLsolve", "SpecialFunctions", "Distributions", "Roots", "ArbNumerics", "StatsBase", "LsqFit", "PoissonRandom", "SparseArrays", "Distributed", "CSV", "SharedArrays", "JLD2","DataFrames","GZip","Parsers","OrderedCollections","FastaIO","RCall","ClusterManagers"])'

ADD scripts/precompile_mktest.jl /root/precompile_mktest.jl
ADD scripts/abcmk_cli.jl /root/abcmk_cli.jl

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/Analytical.jl"))'
RUN julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:Analytical, sysimage_path="/root/mktest.so", precompile_execution_file="/root/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /root/ABCreg
RUN cd /root/ABCreg/src && make

# ENV PATH /root/julia-1.5.1/bin:$PATH
ENV PATH /root/ABCreg/src:$PATH

WORKDIR /analysis

CMD ["julia","--sysimage","/root/mktest.so"]

	