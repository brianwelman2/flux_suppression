#!/bin/bash

cd /home/welman/masters/projects/flux_suppression

# Make a KAT-7 measurement set
rm -r ms/kat7.ms
simms \
--type ascii \
--coord-sys itrf \
--tel kat-7 \
--name kat7.ms \
--outdir ms \
--ra 11h49m15s \
--dec -30d14m50s \
--synthesis-time 2 \
--dtime 10 \
--nchan 1 \
--freq0 1GHz \
--dfreq 1MHz \
--pol XX \
--nolog

# Make a MeerKAT measurement set
rm -r ms/meerkat.ms
simms \
--type ascii \
--coord-sys itrf \
--tel meerkat \
--name meerkat.ms \
--outdir ms \
--ra 11h50m15s \
--dec -30d27m43s \
--synthesis-time 2 \
--dtime 10 \
--nchan 1 \
--freq0 1GHz \
--dfreq 1MHz \
--pol XX \
--nolog

# Make a VLA-B measurement set
rm -r ms/vlab.ms
simms \
--type ascii \
--coord-sys itrf \
--tel vla-b \
--name vlab.ms \
--outdir ms \
--ra 11h49m39s \
--dec 30d08m35s \
--synthesis-time 2 \
--dtime 10 \
--nchan 1 \
--freq0 1GHz \
--dfreq 1MHz \
--pol XX \
--nolog
