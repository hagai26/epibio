#!/bin/tcsh

foreach n (`seq 1 1 32`)
	sbatch --output=tcga_run_out/tcga_${n}.out run_tcga.csh $n
end
