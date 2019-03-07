#!/bin/bash

#SBATCH -J r-oce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -A iPlant-Collabs

module load tacc-singularity

IMG=/work/05066/imicrobe/singularity/c-microbial-map-0.0.1.img
CENTROIDS=/work/05066/imicrobe/iplantc.org/data/scope/centroids/centroids.db
BLAST=/work/05066/imicrobe/iplantc.org/data/scope/blast/blast

set -u

singularity exec $IMG blast2cmap.py "$@" -c "$CENTROIDS" -b "$BLAST" -o cmap-out

echo "Done."
echo "Comments to Ken Youens-Clark <kyclark@email.arizona.edu>"
