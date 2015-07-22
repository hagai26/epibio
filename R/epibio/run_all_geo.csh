#!/bin/tcsh

mkdir -p output/geo_run_out
foreach n (`seq 1 1 300`)
	sbatch --output=output/geo_run_out/geo_${n}.out run_geo.csh $n
end
