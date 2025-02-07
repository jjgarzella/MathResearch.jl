#!/bin/bash

if [ ! -d "ToricControlledReduction" ]; then
  git clone https://github.com/edgarcosta/ToricControlledReduction.git
fi

mkdir data

# sudo apt-get install sagemath

cd controlledreduction

./configure

make
