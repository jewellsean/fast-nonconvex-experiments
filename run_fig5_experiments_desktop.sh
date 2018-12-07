#!/usr/bin/env bash

cd src
for measure in Cor VP VR; do
	for dataset in 2; do
	    for cell_i in 12; do
            Rscript fig5_experiments.R ${measure} ${dataset} ${cell_i} constrained ../configs/configs_fig5.yml
	    done
	done
done