#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J cmaptest

set -u

./run.sh -q "$WORK/data/scope/c-microbial-map/Prochlorococcus_example_query.fasta"
