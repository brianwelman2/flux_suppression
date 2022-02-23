#!/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression

# Image KAT7 data - CLEAN
wsclean -j 8 -size 50 50 -scale 30asec -mgain 0.2 -threshold 0.002 -data-column CLEAN_SRC1_DATA -name images/kat7/src-1/kat7-true-src-1 -niter 1000000 -nmiter 10 -weight uniform ms/kat7.ms
wsclean -j 8 -size 100 100 -scale 30asec -mgain 0.6 -threshold 0.002 -data-column CLEAN_SRC2_DATA -name images/kat7/src-2/kat7-true-src-2 -niter 1000000 -nmiter 10 -weight uniform ms/kat7.ms

Image VLAB data - CLEAN
wsclean -j 8 -size 500 500 -scale 0.5asec -mgain 0.8 -threshold 0.002 -data-column CLEAN_SRC1_DATA -name images/vlab/src-1/vlab-true-src-1 -niter 10000 -weight uniform ms/vlab.ms
wsclean -j 8 -size 500 500 -scale 3asec -mgain 0.8 -threshold 0.002 -data-column CLEAN_SRC2_DATA -name images/vlab/src-2/vlab-true-src-2 -niter 1000000 -nmiter 10 -weight uniform ms/vlab.ms
wsclean -j 8 -size 500 500 -scale 8asec -mgain 0.8 -threshold 0.002 -data-column CLEAN_FP100_SRC100_DATA -name images/vlab/src-100/vlab-true-src-100 -niter 1000000 -nmiter 10 -weight uniform ms/vlab.ms

# Image MeerKAT data - CLEAN
wsclean -j 8 -size 500 500 -scale 0.5asec -mgain 0.2 -threshold 0.002 -data-column CLEAN_SRC1_DATA -name images/meerkat/src-1/meerkat-true-src-1 -niter 1000000 -nmiter 10 -weight uniform ms/meerkat.ms
wsclean -j 8 -size 500 500 -scale 3asec -mgain 0.2 -threshold 0.002 -data-column CLEAN_SRC2_DATA -name images/meerkat/src-2/meerkat-true-src-2 -niter 1000000 -nmiter 10 -weight uniform ms/meerkat.ms
wsclean -j 8 -size 500 500 -scale 8asec -mgain 0.2 -threshold 0.002 -data-column CLEAN_FP100_SRC100_DATA -name images/meerkat/src-100/meerkat-true-src-100 -niter 10000000 -nmiter 10 -weight uniform ms/meerkat.ms