#!/bin/bash
set -e

image="singularity/abcmk.sif"
def="singularity/abcmk.def"

sudo singularity build -s ${image} ${def}
