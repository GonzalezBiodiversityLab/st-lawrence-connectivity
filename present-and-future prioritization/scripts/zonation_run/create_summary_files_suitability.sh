#!/bin/bash
base_circuitscape="/lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/Circuitscape/outputs/"
base_ststim="/lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/STSim/libraries/final/"

#scenarios=(NCNC NC45 NC85 BAUNC BAU45 BAU85 CONNC CON45 CON85)
scenarios=(NC_NC NC_45 NC_85 BAU_NC BAU_45 BAU_85 CON_NC CON_45 CON_85)
species=(BLBR MAAM PLCI RANA URAM)

sc=$2

for timestep in `seq 2010 20 2110`; do
		if [ "$sc" -gt 5 ]; 
		then
			base2_ststim="/BTSL_stconnect.ssim.output/Scenario-166/stconnect_HSOutputHabitatSuitability"
			base3_ststim="/BTSL_stconnect.ssim.output/Scenario-166/stconnect_HSOutputHabitatPatch"
		else
			base2_ststim="/BTSL_stconnect.ssim.output/Scenario-164/stconnect_HSOutputHabitatSuitability"
			base3_ststim="/BTSL_stconnect.ssim.output/Scenario-164/stconnect_HSOutputHabitatPatch"
		fi
		gdal_calc.py --calc="A" -A "${base_ststim}${scenarios[sc]}${base2_ststim}/HabitatSuitability.${species[$1]}.it1.ts${timestep}.tif" --co="COMPRESS=LZW" --overwrite --format=GTiff --outfile="/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.${scenarios[sc]}.${species[$1]}.ts${timestep}.tif"
		for i in `seq 2 1 40`; do
				gdal_calc.py --calc="A+B" -A "/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.${scenarios[sc]}.${species[$1]}.ts${timestep}.tif" -B "${base_ststim}${scenarios[sc]}${base2_ststim}/HabitatSuitability.${species[$1]}.it${i}.ts${timestep}.tif" --co="COMPRESS=LZW" --overwrite --format=GTiff --outfile="/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.${scenarios[sc]}.${species[$1]}.ts${timestep}.tif"
		done
done

