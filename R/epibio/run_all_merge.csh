#!/bin/tcsh

mkdir -p output/merge_run_out
foreach n (`seq 1 1 200`)
	sbatch --output=output/merge_run_out/geo_${n}.out run_merge.csh $n
end
