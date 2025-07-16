#!/bin/bash
echo $1
out=/lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/outputs_btls/${1}.ABF_ME.rank.compressed.tif
if test -f "$out"; then
	echo "$out exists"
else
	/home/glaroc/zonation/build/zig4/zig4 -r ../inputs/settings.dat ${1} ../outputs_btsl/${1}.txt 0.0 0 1.0 0 0 --grid-output-formats compressed-tif
fi