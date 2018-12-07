#!/usr/bin/env bash

echo "Creating Figures "
cd src
Rscript fig5_bottom_create.R ../configs/configs_fig5.yml
Rscript fig5_top_create.R ../configs/configs_fig5_bottom.yml
