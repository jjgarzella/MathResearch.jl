#!/bin/bash

# source ~/mambaforge/etc/profile.d/conda.sh
# conda activate toriccr

./controlledreduction/build/examples/readfile $1 data/$(basename $1).out
