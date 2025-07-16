#!/bin/bash
species=(BLBR MAAM PLCI RANA URAM)

for sp in `seq 0 4`; do
	gdal_calc.py --calc="A-B" -A "/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.NC_85.${species[sp]}.ts2010.tif" -B "/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.NC_85.${species[sp]}.ts2110.tif" --co="COMPRESS=LZW" --overwrite --format=GTiff --outfile="/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.NC85.${species[sp]}.2110-2010.tif"


	gdal_calc.py --calc="abs(A-B)" -A "/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.BAU_NC.${species[sp]}.ts2010.tif" -B "/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.BAU_NC.${species[sp]}.ts2110.tif" --co="COMPRESS=LZW" --overwrite --format=GTiff --outfile="/home/glaroc/projects/rrg-gonzalez/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/stsim_summary_maps/HabitatSuitability.BAUNC.${species[sp]}.2110-2010.tif"
done