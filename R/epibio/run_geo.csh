#!/bin/csh

if ($#argv != 1) then
	echo "Usage: $0 <# GEO>"
		exit 0
endif
source /cs/icore/joshua.moss/dor/hagaic/epibio/hagaic_env.sh
#SBATCH --mem=128000
#SBATCH --time=72:00:00
srun Rscript organize_geo_jm.R ${1}


