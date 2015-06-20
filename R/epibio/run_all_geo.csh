#!/bin/tcsh

foreach n (`seq 1 1 100`)
	sbatch --output=output/geo_run_out/geo_${n}.out run_geo.csh $n
end
