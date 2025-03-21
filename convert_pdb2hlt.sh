#!/bin/bash

original_pdb="$1"

# get the directory and name of the original pdb
original_dir=$(dirname "$original_pdb")
name=$(basename "$original_pdb")
name="${name%.*}"

# set the path and name of the hlt pdb
hlt_dir="${original_dir}/HLT"
mkdir -p "$hlt_dir"

hlt_pdb="${hlt_dir}/${name}.pdb"

# convert pdb to hlt pdb
poetry run python /home/scripts/util/chothia2HLT.py "$original_pdb" --heavy H --light L --target A --output "$hlt_pdb"

echo "Converting finished: $hlt_pdb !"