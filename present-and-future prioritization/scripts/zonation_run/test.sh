#!/bin/bash
scenarios=(NCNC NC45 NC85 BAUNC BAU45 BAU85 CONNC CON45 CON85)
sce=($(seq 0 8))
it=($(seq 2010 20 2110))
sce_ts=()
for sc in ${scenarios[@]}; do
	for timestep in `seq 2010 20 2110`; do
		echo $sc
		echo $timestep
	done
done
