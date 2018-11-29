#!/usr/bin/env bash

cd src
for measure in Cor VP VR; do
	for dataset in 2; do
	    for cell_i in 12; do
            Rscript fig6_experiments.R ${measure} ${dataset} ${cell_i} constrained ../configs/configs_fig6.yml
	    done
	done
done