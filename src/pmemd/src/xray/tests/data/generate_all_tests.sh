#!/bin/bash

# Requires proper conda env to be active, e.g:
# source ~/miniconda3/bin/activate
# Install reference cctbx package
# conda install -c conda-forge cctbx=2021.7
conda run python gen_scaling_atomic_test.py
conda run python get_ml_alpha_beta_in_zones_test.py
conda run python get_ml_alpha_beta_individual_test.py
conda run python get_scaling_log_binning_test.py