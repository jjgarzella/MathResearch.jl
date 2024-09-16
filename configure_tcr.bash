#!/bin/bash

if [ ! -d "ToricControlledReduction" ]; then
  git clone https://github.com/edgarcosta/ToricControlledReduction.git
fi

mkdir data

# some conda config possibly necessary
# not really sure of the best way to install conda

# the C_INCLUDE_PATH for sure and maybe LD_LIBRARY_PATH need to be set for it to work

# conda create --name toriccr ntl
# conda env config vars set C_INCLUDE_PATH=/u/jg4411/mambaforge/envs/toriccr/include:${C_INCLUDE_PATH} -n toriccr
# conda env config vars set LD_LIBRARY_PATH=/u/jg4411/mambaforge/envs/toriccr/lib:${LD_LIBRARY_PATH} -n toriccr

# conda activate toriccr

cd ToricControlledReduction

./configure # --with-ntl=~/mambaforge/envs/toriccr
# ./configure --with-ntl=/u/jg4411/mambaforge/envs/toriccr

# make check
make examples

