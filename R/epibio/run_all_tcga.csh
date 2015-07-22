#!/bin/tcsh

mkdir -p output/tcga_run_out
foreach n (`seq 1 1 50`)
	sbatch --output=output/tcga_run_out/tcga_${n}.out run_tcga.csh $n
end
