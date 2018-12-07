#!/usr/bin/env bash

cd src
for T in 100 215 464 1000 2154 4641 10000 21544 46415 100000; do
	for THETA in 0.1 0.01 0.001; do
        Rscript figS2-3_experiments.R ${T} ${THETA} ../configs/configs_figS2-3.yml
	done
done
