#!/bin/bash
base_circuitscape="/lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/Circuitscape/outputs/"
base_ststim="/lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/STSim/libraries/final/"

scenarios=(NCNC NC45 NC85 BAUNC BAU45 BAU85 CONNC CON45 CON85)
scenarios2=(NC_NC NC_45 NC_85 BAU_NC BAU_45 BAU_85 CON_NC CON_45 CON_85)
species=(BLBR MAAM PLCI RANA URAM)

for sc in {0..8}; do
	for timestep in `seq 2010 20 2110`; do
		for i in {1..40}; do
			if [ "$sc" -gt 5 ]; 
			then
				base2_ststim="/BTSL_stconnect.ssim.output/Scenario-166/stconnect_HSOutputHabitatSuitability"
				base3_ststim="/BTSL_stconnect.ssim.output/Scenario-166/stconnect_HSOutputHabitatPatch"
			else
				base2_ststim="/BTSL_stconnect.ssim.output/Scenario-164/stconnect_HSOutputHabitatSuitability"
				base3_ststim="/BTSL_stconnect.ssim.output/Scenario-164/stconnect_HSOutputHabitatPatch"
			fi
			fn="/lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/inputs/spp/${scenarios[${sc}]}_${timestep}_it${i}.spp"
			touch ${fn}
			for s in "${species[@]}"; do
				printf "1 0.04 1 1 0.25 ${base_ststim}${scenarios2[${sc}]}${base2_ststim}/HabitatSuitability.${s}.it${i}.ts${timestep}.tif\n" >> ${fn}
				printf "1 0.04 1 1 0.25 ${base_circuitscape}${scenarios[${sc}]}/curmap_${scenarios[${sc}]}_Resistance.${s}.it${i}.ts${timestep}.tif\n" >> ${fn}
				printf "1 0.04 1 1 0.25 ${base_ststim}${scenarios2[${sc}]}${base3_ststim}/HabitatPatch.${s}.it${i}.ts${timestep}.tif\n" >> ${fn}
			done
		done
	done
done
