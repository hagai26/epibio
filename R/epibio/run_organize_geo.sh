#!/bin/sh
#SBATCH --mem=64000
#SBATCH --output=organize_geo_sbatch.out
#SBATCH --time=72:00:00
srun Rscript organize_geo.R
