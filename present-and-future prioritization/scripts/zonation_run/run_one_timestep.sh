#!/bin/bash

scenarios=(NCNC NC45 NC85 BAUNC BAU45 BAU85 CONNC CON45 CON85)
sce=($(seq 0 8))
it=($(seq 2010 20 2110))
sce_ts=()
for sc in ${scenarios[@]}; do
	for timestep in `seq 2010 20 2110`; do
			sce_ts+=("${sc}_${timestep}")
	done
done
name=${sce_ts[$1]}'*.spp'
echo $name
find /lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/inputs/spp/ -name $name | parallel -j 40 /lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/run/run_one_iteration.sh {1} 

