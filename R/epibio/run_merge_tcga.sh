#!/bin/sh
#SBATCH --mem=64000
#SBATCH --output=merge_tcga.out 
#SBATCH --time=72:00:00
srun Rscript merge_beta_values_tcga.R
