#!/bin/tcsh

foreach n (`seq 1 1 157`)
	sbatch --output=merge_run_out/geo_${n}.out run_merge.csh $n
end
