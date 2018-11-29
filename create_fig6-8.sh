#!/usr/bin/env bash

echo "Creating Figures "
cd src
Rscript fig6_create.R ../configs/configs_fig6.yml
Rscript fig8_create.R ../configs/configs_fig8.yml
