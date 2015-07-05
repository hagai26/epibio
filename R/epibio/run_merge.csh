#!/bin/csh
#SBATCH --mem=128000
#SBATCH --time=72:00:00

if ($#argv != 1) then
	echo "Usage: $0 <# GEO>"
		exit 0
endif
source /cs/icore/joshua.moss/dor/hagaic/epibio/hagaic_env_3.2.1.sh
srun Rscript merge_beta_values.R ${1}


