#!/bin/bash
set -e

image="singularity/abcmk.simg"
def="singularity/abcmk.def"

sudo singularity build -s ${image} ${def}
sudo chown ${USER}:${image}
