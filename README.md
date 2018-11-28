# fast-nonconvex-experiments

This repository contains code to reproduce all experiments in Jewell, Hocking, Fearnhead, and Witten (2018). [Fast Nonconvex Deconvolution of Calcium Imaging Data](https://arxiv.org/abs/1802.07380). 

# Requires:

Experiments require the following two packages 

- [FastLZeroSpikeInference](https://github.com/jewellsean/FastLZeroSpikeInference)
- [OASIS](https://github.com/j-friedrich/OASIS)

where the R implementation of FastLZeroSpikeInference should be setup following [these instructions](https://github.com/jewellsean/FastLZeroSpikeInference) and OASIS should be setup inside a conda environment named ```aibs_py35```. To setup this environment with OASIS create the following environment file named ```aibs_py35.yml```

```
name: aibs_py35
dependencies:
  - python=3.5  
  - matplotlib
  - numpy
  - scipy
  - cython=0.23.4
  - pandas
```

Then run 
```
conda env create -f aibs_py35.yml
```
to setup the conda environment. 

To install OASIS into this directory, run

```
cd OASIS

source activate aibs_py35.yml

python3 setup.py build_ext --inplace
python3 setup.py install
python3 setup.py clean --all
```

