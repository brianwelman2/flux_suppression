#!/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression

# # Make a KAT-7 gains for each skymodel (1 and 2)
rm gains/true/kat7/*.npy
kalcal configs/kat7/simulate_kat7_gains.yml

# # Make a MeerKAT gains for each skymodel (1, 2 and 100)
rm gains/true/meerkat/*.npy
kalcal configs/meerkat/simulate_meerkat_gains.yml

# Make a VLA-B gains for each skymodel (1, 2 and 100)
rm gains/true/vlab/*.npy
kalcal configs/vlab/simulate_vlab_gains.yml