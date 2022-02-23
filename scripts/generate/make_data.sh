#!/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression

# Make a KAT-7 model visibilities, true visibilities and noisy visibilities
kalcal configs/kat7/simulate_kat7_data.yml

# Make a MeerKAT model visibilities, true visibilities and noisy visibilities
kalcal configs/meerkat/simulate_meerkat_data.yml

# Make a VLA-B model visibilities, true visibilities and noisy visibilities
kalcal configs/vlab/simulate_vlab_data.yml