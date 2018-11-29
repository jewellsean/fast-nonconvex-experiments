#!/bin/bash
source ~/.bashrc
source activate aibs_py35

echo $(which python)
python utils/run_oasis.py "$@"