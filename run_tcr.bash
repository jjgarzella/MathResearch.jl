#!/bin/bash

# source ~/mambaforge/etc/profile.d/conda.sh
# conda activate toriccr

./ToricControlledReduction/build/examples/readfile.exe $1 data/$(basename $1).out
