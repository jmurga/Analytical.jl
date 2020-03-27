# Analytical

[![Build Status](https://travis-ci.com/jmurga/Analytical.jl.svg?branch=master)](https://travis-ci.com/jmurga/Analytical.jl)

Extended ABC-MK calculations accounting for background selection and weak adaptation based on Julia language.

Included in this repository are the scripts used to simulate and infer parameters from ABC software.

Debian imaging including *julia* and all packages need to run the estimations. It also includes *jupyterlab* and several data science python packages. To run the images (jupyter notebook will start on localhost:8888):

```bash
# Pull the image
docker pull jmurga/mktest
# Run docker bash interactive session
docker run -i -t -p 8888:8888 -v ${HOME}/<anyData>:/data/<anyData> -e HOSTID=$(id -u) jmurga/mktest
# Run only jupyter notebook from docker image
docker run -i -t -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter notebook --ip='*' --port=8888 --no-browser"
```
